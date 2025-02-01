// tools::package_native_routine_registration_skeleton("fbati")
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void condGeneFBATControl_backupTrait(void *, void *);
extern void condGeneFBATControl_centerTrait(void *, void *, void *);
extern void condGeneFBATControl_contsImc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void condGeneFBATControl_contsUimc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void condGeneFBATControl_dUdBc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void condGeneFBATControl_estEq(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void condGeneFBATControl_estEqNuis(void *, void *, void *, void *, void *);
extern void condGeneFBATControl_estEqNuisUpdate(void *, void *, void *);
extern void condGeneFBATControl_estEqNuisUpdate2(void *, void *, void *);
extern void condGeneFBATControl_free(void *);
extern void condGeneFBATControl_imc(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void condGeneFBATControl_linkTrait(void *, void *, void *, void *);
extern void condGeneFBATControl_load(void *, void *);
extern void condGeneFBATControl_numFam(void *, void *);
extern void condGeneFBATControl_numInfFam(void *, void *);
extern void condGeneFBATControl_pids(void *, void *);
extern void condGeneFBATControl_print(void *);
extern void condGeneFBATControl_proportionInformative(void *, void *);
extern void condGeneFBATControl_removeUnphased(void *);
extern void condGeneFBATControl_restoreTrait(void *, void *);
extern void condGeneFBATControl_robustStat(void *, void *, void *, void *, void *, void *);
extern void condGeneFBATControl_uimc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void condGeneFBATControl_varContsMean(void *, void *, void *, void *);
extern void condGeneFBATControl_varContsModel(void *, void *, void *, void *);
extern void condGeneFBATControl_varExplConts(void *, void *, void *, void *);
//extern void cpp_gesim_clear();
extern void cpp_gesim_clear(void);
//extern void cpp_gesim_draw();
extern void cpp_gesim_draw(void);
//extern void cpp_gesim_print();
extern void cpp_gesim_print(void);
extern void cpp_gesim_set(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
//extern void cpp_gped_clear();
extern void cpp_gped_clear(void);
extern void cpp_gped_estEq(void *, void *, void *);
extern void cpp_gped_numCovariates(void *);
extern void cpp_gped_print(void *);
extern void cpp_gped_set(void *, void *, void *, void *, void *, void *, void *, void *);
extern void cpp_gped_set_C(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void cpp_gped_setStrategy(void *);
extern void cpp_gped_statCompute(void *, void *, void *);
extern void cpp_gped_statCompute_A(void *, void *, void *);
extern void cpp_gped_statPrecompute(void *, void *, void *);
extern void cpp_gped_statPrecompute_A(void *, void *, void *);
//extern void cpp_mmatrix_debug();
extern void cpp_mmatrix_debug(void);
//extern void cpp_rn_attach();
extern void cpp_rn_attach(void);
//extern void cpp_rn_debug();
extern void cpp_rn_debug(void);
//extern void cpp_rn_detach();
extern void cpp_rn_detach(void);
extern void cpp_rn_setNormalSigma(void *, void *);
extern void ddataComputeGroupG(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

//extern void eREXP_fbatme(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
//extern void eREXP_fbatmeev(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
//extern void eREXP_joint(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

// 2020-07-07: LTO error fixes -- too many args?
extern void eREXP_fbatme(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void eREXP_fbatmeev(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void eREXP_joint(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);




//extern void fbati_cpp(void *, void *, void *, void *, void *, void *); // 2025-01-05
extern void nuclify_cpp(void *, void *, void *, void *, void *);
extern void pG_group_dehash(void *, void *);
extern void strataReduce_cpp(void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"condGeneFBATControl_backupTrait",           (DL_FUNC) &condGeneFBATControl_backupTrait,            2},
    {"condGeneFBATControl_centerTrait",           (DL_FUNC) &condGeneFBATControl_centerTrait,            3},
    {"condGeneFBATControl_contsImc",              (DL_FUNC) &condGeneFBATControl_contsImc,              10},
    {"condGeneFBATControl_contsUimc",             (DL_FUNC) &condGeneFBATControl_contsUimc,             11},
    {"condGeneFBATControl_dUdBc",                 (DL_FUNC) &condGeneFBATControl_dUdBc,                 13},
    {"condGeneFBATControl_estEq",                 (DL_FUNC) &condGeneFBATControl_estEq,                  9},
    {"condGeneFBATControl_estEqNuis",             (DL_FUNC) &condGeneFBATControl_estEqNuis,              5},
    {"condGeneFBATControl_estEqNuisUpdate",       (DL_FUNC) &condGeneFBATControl_estEqNuisUpdate,        3},
    {"condGeneFBATControl_estEqNuisUpdate2",      (DL_FUNC) &condGeneFBATControl_estEqNuisUpdate2,       3},
    {"condGeneFBATControl_free",                  (DL_FUNC) &condGeneFBATControl_free,                   1},
    {"condGeneFBATControl_imc",                   (DL_FUNC) &condGeneFBATControl_imc,                    9},
    {"condGeneFBATControl_linkTrait",             (DL_FUNC) &condGeneFBATControl_linkTrait,              4},
    {"condGeneFBATControl_load",                  (DL_FUNC) &condGeneFBATControl_load,                   2},
    {"condGeneFBATControl_numFam",                (DL_FUNC) &condGeneFBATControl_numFam,                 2},
    {"condGeneFBATControl_numInfFam",             (DL_FUNC) &condGeneFBATControl_numInfFam,              2},
    {"condGeneFBATControl_pids",                  (DL_FUNC) &condGeneFBATControl_pids,                   2},
    {"condGeneFBATControl_print",                 (DL_FUNC) &condGeneFBATControl_print,                  1},
    {"condGeneFBATControl_proportionInformative", (DL_FUNC) &condGeneFBATControl_proportionInformative,  2},
    {"condGeneFBATControl_removeUnphased",        (DL_FUNC) &condGeneFBATControl_removeUnphased,         1},
    {"condGeneFBATControl_restoreTrait",          (DL_FUNC) &condGeneFBATControl_restoreTrait,           2},
    {"condGeneFBATControl_robustStat",            (DL_FUNC) &condGeneFBATControl_robustStat,             6},
    {"condGeneFBATControl_uimc",                  (DL_FUNC) &condGeneFBATControl_uimc,                  12},
    {"condGeneFBATControl_varContsMean",          (DL_FUNC) &condGeneFBATControl_varContsMean,           4},
    {"condGeneFBATControl_varContsModel",         (DL_FUNC) &condGeneFBATControl_varContsModel,          4},
    {"condGeneFBATControl_varExplConts",          (DL_FUNC) &condGeneFBATControl_varExplConts,           4},
    {"cpp_gesim_clear",                           (DL_FUNC) &cpp_gesim_clear,                            0},
    {"cpp_gesim_draw",                            (DL_FUNC) &cpp_gesim_draw,                             0},
    {"cpp_gesim_print",                           (DL_FUNC) &cpp_gesim_print,                            0},
    {"cpp_gesim_set",                             (DL_FUNC) &cpp_gesim_set,                             21},
    {"cpp_gped_clear",                            (DL_FUNC) &cpp_gped_clear,                             0},
    {"cpp_gped_estEq",                            (DL_FUNC) &cpp_gped_estEq,                             3},
    {"cpp_gped_numCovariates",                    (DL_FUNC) &cpp_gped_numCovariates,                     1},
    {"cpp_gped_print",                            (DL_FUNC) &cpp_gped_print,                             1},
    {"cpp_gped_set",                              (DL_FUNC) &cpp_gped_set,                               8},
    {"cpp_gped_set_C",                            (DL_FUNC) &cpp_gped_set_C,                            10},
    {"cpp_gped_setStrategy",                      (DL_FUNC) &cpp_gped_setStrategy,                       1},
    {"cpp_gped_statCompute",                      (DL_FUNC) &cpp_gped_statCompute,                       3},
    {"cpp_gped_statCompute_A",                    (DL_FUNC) &cpp_gped_statCompute_A,                     3},
    {"cpp_gped_statPrecompute",                   (DL_FUNC) &cpp_gped_statPrecompute,                    3},
    {"cpp_gped_statPrecompute_A",                 (DL_FUNC) &cpp_gped_statPrecompute_A,                  3},
    {"cpp_mmatrix_debug",                         (DL_FUNC) &cpp_mmatrix_debug,                          0},
    {"cpp_rn_attach",                             (DL_FUNC) &cpp_rn_attach,                              0},
    {"cpp_rn_debug",                              (DL_FUNC) &cpp_rn_debug,                               0},
    {"cpp_rn_detach",                             (DL_FUNC) &cpp_rn_detach,                              0},
    {"cpp_rn_setNormalSigma",                     (DL_FUNC) &cpp_rn_setNormalSigma,                      2},
    {"dataComputeGroupG",                         (DL_FUNC) &ddataComputeGroupG,                         11},
    //{"eREXP_fbatme",                              (DL_FUNC) &eREXP_fbatme,                              11},
    //{"eREXP_fbatmeev",                            (DL_FUNC) &eREXP_fbatmeev,                            11},
    //{"eREXP_joint",                               (DL_FUNC) &eREXP_joint,                               14},
    {"eREXP_fbatme",                              (DL_FUNC) &eREXP_fbatme,                              10},
    {"eREXP_fbatmeev",                            (DL_FUNC) &eREXP_fbatmeev,                            10},
    {"eREXP_joint",                               (DL_FUNC) &eREXP_joint,                               13},
    //{"fbati_cpp",                                 (DL_FUNC) &fbati_cpp,                                   6},  // commented out 2025-01-05
    {"nuclify_cpp",                               (DL_FUNC) &nuclify_cpp,                                5},
    {"pG_group_dehash",                           (DL_FUNC) &pG_group_dehash,                            2},
    {"strataReduce_cpp",                          (DL_FUNC) &strataReduce_cpp,                           8},
    {NULL, NULL, 0}
};

void R_init_fbati(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
