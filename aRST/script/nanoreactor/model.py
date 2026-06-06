"""
Q(l+1) = Q(l) + tau * (F - D)
Q_s: population of specie
abundance: normalized population ->dict
F_s: tot. formation flux of specie
D_s: tot. consume flux of specie
v_s:
f_i: forward path flux         for A+A->P: fi = k[C_A]^2/2
Counting all reaction path, sum over all consumption flux for reactant side and production flux for product side
F_s = v_s * generation flux
D_s = c_s * consumption flux
tau: consumption cutoff, species s can only be consumed no more than 50% in this layer

critieia for suriving spiese for next layer: 1. relative population cutoff; 2. top-N population


"""
import numpy as np
from collections import defaultdict, Counter
import math

R_KCAL = 0.0019872041  # kcal mol-1 K-1


def safe_exp(x, lower=-1e3, upper=1e3):
    return np.exp(max(lower, min(upper, x)))


def normalize_population(population):
    total = sum(v for v in population.values() if v > 0.0)
    if total <= 0.0:
        return {}

    return {
        sid: q / total
        for sid, q in population.items()
        if q > 0.0
    }


def safe_TSbarrier(dG_ts):
    # set the lowest TS barrier for barrierless rxn = 1kcal/mol
    if not np.isinf(dG_ts):
        return max(dG_ts, 1)


def get_C(species_counter, abundance):
    # A + B - > X + Y
    # C_products = C_X * C_Y

    val = 1.0

    for sid, coeff in species_counter.items():
        c = abundance.get(sid, 0.0)

        if c <= 0.0:
            return 0.0

        val *= c ** coeff

    return val


def symmetry_factor(species_counter):
    """
    Prevent overcounting identical reactants.

    A + B      sf = 1
    A + A      sf = 2
    A + A + A  sf = 6
    A + A + B  sf = 2
    """

    fac = 1.0

    for coeff in species_counter.values():
        fac *= math.factorial(coeff)

    return fac


def path_flux(info, abundance, temp=298.15, evalmode="g",elim=50):
    """
    :param info: rxn summarize->dict
    :param abundance: of each structure ->dict
    :param temp: of rxn
    :return: r:flux of this rxn, max(r_net, 0) for the forward exploration

    r_net = k_f C_products - k_r C_reactants
    """

    reactants = info.get("reactants", [])
    products = info.get("products", [])
    badP = ["unknow","-","?",None]

    if not reactants:
        return 0.0

    if any(x in badP for x in products):
        return 0.0


    if info.get("bond_change") is False:
        return 0.0

    if evalmode == "g":
        dG_ts = info.get("dG_ts", np.inf)
        G_rxn = info.get("G_rxn", np.inf)
    else:
        dG_ts = info.get("dE_ts", np.inf)
        G_rxn = info.get("E_rxn", np.inf)

    if any(x == np.inf or x>elim for x in [dG_ts, G_rxn]):
        return 0.0

    RT = R_KCAL * temp
    # forward barrier
    dG_ts_f = safe_TSbarrier(dG_ts)

    # reverse barrier
    dG_ts_r = safe_TSbarrier(dG_ts_f - G_rxn)

    k_f = safe_exp(-dG_ts_f / RT)
    k_r = safe_exp(-dG_ts_r / RT)

    r_counter = Counter(reactants)
    p_counter = Counter(products)

    C_R = get_C(r_counter, abundance)
    C_P = get_C(p_counter, abundance)

    if C_R <= 0.0:
        return 0.0

    sf_r = symmetry_factor(r_counter)
    sf_p = symmetry_factor(p_counter)

    forward_flux = k_f * C_R / sf_r
    reverse_flux = k_r * C_P / sf_p

    net_flux = forward_flux - reverse_flux

    return max(net_flux, 0.0)


def consume_cutoff(population, consumption_flux, max_consumed_fraction=0.50):
    """
    for each spieces, can maxmal 50% consumed in this layer
    tau*D_s/Q_s < 0.5
    tau <= 0.5 * Q_s / D_s
    Q: old population from last layer
    D: consumed flux of this layer
    """
    tau_candidates = []
    for sid, D in consumption_flux.items():
        Q = population.get(sid, 0.0)

        if Q > 0.0 and D > 0.0:
            tau_candidates.append(max_consumed_fraction * Q / D)

    if not tau_candidates:
        return 1.0

    return min(tau_candidates)


def accumulate_layer_flux(print_info, abundance, temp=298.15, evalmode="g",elim=50):
    production_flux = defaultdict(float)
    consumption_flux = defaultdict(float)

    path_fluxes = []
    for info in print_info:
        flux = path_flux(
            info,
            abundance=abundance,
            temp=temp,
            evalmode=evalmode,
            elim=elim
        )

        if flux <= 0.0:
            # forward exploration
            continue

        reactants = info.get("reactants", [])
        products = info.get("products", [])

        path_fluxes.append((info, flux))

        # consume
        for r, coeff in Counter(reactants).items():
            consumption_flux[r] += coeff * flux

        # formation
        for p, coeff in Counter(products).items():
            if p not in {"?", "unknow", None, "-"}:
                production_flux[p] += coeff * flux

    return production_flux, consumption_flux, path_fluxes


def update_population(population_old, production_flux, consumption_flux,
                      keep_top=100, rel_cut=1e-5, max_consumed_fraction=0.50):
    tau = consume_cutoff(population_old, consumption_flux, max_consumed_fraction)
    population_raw = defaultdict(float)
    all_species = (
            set(population_old)
            | set(production_flux)
            | set(consumption_flux)
    )
    for sid in all_species:
        old_q = population_old.get(sid, 0.0)
        F = production_flux.get(sid, 0.0)
        D = consumption_flux.get(sid, 0.0)

        q_new = old_q + tau * (F - D)
        population_raw[sid] = max(q_new, 0.0)

    # remove 0 species
    population_positive = {
        sid: q
        for sid, q in population_raw.items()
        if q > 0.0
    }

    if not population_positive:
        return {}, {}, tau, {}

    # relative cutoffc
    max_q = max(population_positive.values())
    population_filtered = {
        sid: q
        for sid, q in population_positive.items()
        if q / max_q >= rel_cut
    }

    # top-N cutoff
    ranked = sorted(
        population_filtered.items(),
        key=lambda x: x[1],
        reverse=True,
    )

    ranked = ranked[:keep_top]

    population_new = dict(ranked)
    abundance_new = normalize_population(population_new)
    removed_info = {}

    old_species = set(population_old)
    new_species = set(population_new)
    raw_species = set(population_raw)

    for sid in old_species | raw_species:
        if sid in new_species:
            continue

        q_raw = population_raw.get(sid, 0.0)

        if q_raw <= 0.0:
            reason = "depleted_by_flux"
        elif q_raw / max_q < rel_cut:
            reason = "below_relative_cutoff"
        else:
            reason = "removed_by_topN_pruning"

        removed_info[sid] = {
            "reason": reason,
            "raw_population": q_raw,
            "old_population": population_old.get(sid, 0.0),
            "production_flux": production_flux.get(sid, 0.0),
            "consumption_flux": consumption_flux.get(sid, 0.0),
        }

    return population_new, abundance_new, tau, removed_info
