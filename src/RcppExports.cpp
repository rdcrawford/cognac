// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/cognac.h"
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// CalcAlgnSubMatrix
Rcpp::NumericMatrix CalcAlgnSubMatrix(std::vector< std::string > seqs);
RcppExport SEXP _cognac_CalcAlgnSubMatrix(SEXP seqsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector< std::string > >::type seqs(seqsSEXP);
    rcpp_result_gen = Rcpp::wrap(CalcAlgnSubMatrix(seqs));
    return rcpp_result_gen;
END_RCPP
}
// CalcAlgnPartitionDists
std::list < Rcpp::NumericMatrix > CalcAlgnPartitionDists(std::string msaPath, std::string method, std::vector< int > genePartitions);
RcppExport SEXP _cognac_CalcAlgnPartitionDists(SEXP msaPathSEXP, SEXP methodSEXP, SEXP genePartitionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type msaPath(msaPathSEXP);
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< std::vector< int > >::type genePartitions(genePartitionsSEXP);
    rcpp_result_gen = Rcpp::wrap(CalcAlgnPartitionDists(msaPath, method, genePartitions));
    return rcpp_result_gen;
END_RCPP
}
// ConcatenateAlignments
void ConcatenateAlignments(Rcpp::StringVector& concatAlgn, const Rcpp::StringVector& algn);
RcppExport SEXP _cognac_ConcatenateAlignments(SEXP concatAlgnSEXP, SEXP algnSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector& >::type concatAlgn(concatAlgnSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type algn(algnSEXP);
    ConcatenateAlignments(concatAlgn, algn);
    return R_NilValue;
END_RCPP
}
// CreateAlgnDistMat
Rcpp::NumericMatrix CreateAlgnDistMat(std::string msaPath, std::string method);
RcppExport SEXP _cognac_CreateAlgnDistMat(SEXP msaPathSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type msaPath(msaPathSEXP);
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(CreateAlgnDistMat(msaPath, method));
    return rcpp_result_gen;
END_RCPP
}
// CreateCognacRunData
void CreateCognacRunData(Rcpp::Environment& geneEnv, const std::vector< std::string >& gfPaths, const std::vector< std::string >& faPaths, const std::string& faaPath);
RcppExport SEXP _cognac_CreateCognacRunData(SEXP geneEnvSEXP, SEXP gfPathsSEXP, SEXP faPathsSEXP, SEXP faaPathSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Environment& >::type geneEnv(geneEnvSEXP);
    Rcpp::traits::input_parameter< const std::vector< std::string >& >::type gfPaths(gfPathsSEXP);
    Rcpp::traits::input_parameter< const std::vector< std::string >& >::type faPaths(faPathsSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type faaPath(faaPathSEXP);
    CreateCognacRunData(geneEnv, gfPaths, faPaths, faaPath);
    return R_NilValue;
END_RCPP
}
// CreateCoreGenomeDistMat
Rcpp::NumericMatrix CreateCoreGenomeDistMat(std::string msaPath);
RcppExport SEXP _cognac_CreateCoreGenomeDistMat(SEXP msaPathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type msaPath(msaPathSEXP);
    rcpp_result_gen = Rcpp::wrap(CreateCoreGenomeDistMat(msaPath));
    return rcpp_result_gen;
END_RCPP
}
// DeletePartitions
void DeletePartitions(std::string msaPath, std::vector<int> delStart, std::vector<int> delEnd, std::string outPath);
RcppExport SEXP _cognac_DeletePartitions(SEXP msaPathSEXP, SEXP delStartSEXP, SEXP delEndSEXP, SEXP outPathSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type msaPath(msaPathSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type delStart(delStartSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type delEnd(delEndSEXP);
    Rcpp::traits::input_parameter< std::string >::type outPath(outPathSEXP);
    DeletePartitions(msaPath, delStart, delEnd, outPath);
    return R_NilValue;
END_RCPP
}
// ExtractGenomeNameFromPath
std::string ExtractGenomeNameFromPath(std::string path);
RcppExport SEXP _cognac_ExtractGenomeNameFromPath(SEXP pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type path(pathSEXP);
    rcpp_result_gen = Rcpp::wrap(ExtractGenomeNameFromPath(path));
    return rcpp_result_gen;
END_RCPP
}
// GetGenomeNameWithExt
std::string GetGenomeNameWithExt(std::string path, std::string ext);
RcppExport SEXP _cognac_GetGenomeNameWithExt(SEXP pathSEXP, SEXP extSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type path(pathSEXP);
    Rcpp::traits::input_parameter< std::string >::type ext(extSEXP);
    rcpp_result_gen = Rcpp::wrap(GetGenomeNameWithExt(path, ext));
    return rcpp_result_gen;
END_RCPP
}
// FilterAlgnPositions
void FilterAlgnPositions(std::string msaPath, std::string filterMsaPath, double minGapFrac, int minSubThresh);
RcppExport SEXP _cognac_FilterAlgnPositions(SEXP msaPathSEXP, SEXP filterMsaPathSEXP, SEXP minGapFracSEXP, SEXP minSubThreshSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type msaPath(msaPathSEXP);
    Rcpp::traits::input_parameter< std::string >::type filterMsaPath(filterMsaPathSEXP);
    Rcpp::traits::input_parameter< double >::type minGapFrac(minGapFracSEXP);
    Rcpp::traits::input_parameter< int >::type minSubThresh(minSubThreshSEXP);
    FilterAlgnPositions(msaPath, filterMsaPath, minGapFrac, minSubThresh);
    return R_NilValue;
END_RCPP
}
// FilterPartitionedAlgnPositions
std::vector< int > FilterPartitionedAlgnPositions(std::string msaPath, std::string filterMsaPath, std::vector<int> genePositions, double minGapFrac, int minSubThresh);
RcppExport SEXP _cognac_FilterPartitionedAlgnPositions(SEXP msaPathSEXP, SEXP filterMsaPathSEXP, SEXP genePositionsSEXP, SEXP minGapFracSEXP, SEXP minSubThreshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type msaPath(msaPathSEXP);
    Rcpp::traits::input_parameter< std::string >::type filterMsaPath(filterMsaPathSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type genePositions(genePositionsSEXP);
    Rcpp::traits::input_parameter< double >::type minGapFrac(minGapFracSEXP);
    Rcpp::traits::input_parameter< int >::type minSubThresh(minSubThreshSEXP);
    rcpp_result_gen = Rcpp::wrap(FilterPartitionedAlgnPositions(msaPath, filterMsaPath, genePositions, minGapFrac, minSubThresh));
    return rcpp_result_gen;
END_RCPP
}
// FindIdenticalGenes
Rcpp::List FindIdenticalGenes(const Rcpp::StringVector& genes, const Rcpp::StringVector& geneIds);
RcppExport SEXP _cognac_FindIdenticalGenes(SEXP genesSEXP, SEXP geneIdsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type genes(genesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type geneIds(geneIdsSEXP);
    rcpp_result_gen = Rcpp::wrap(FindIdenticalGenes(genes, geneIds));
    return rcpp_result_gen;
END_RCPP
}
// GetAlgnQualScores
std::vector< double > GetAlgnQualScores(std::string msaPath, std::string method, int stepVal, int windowSize);
RcppExport SEXP _cognac_GetAlgnQualScores(SEXP msaPathSEXP, SEXP methodSEXP, SEXP stepValSEXP, SEXP windowSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type msaPath(msaPathSEXP);
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< int >::type stepVal(stepValSEXP);
    Rcpp::traits::input_parameter< int >::type windowSize(windowSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(GetAlgnQualScores(msaPath, method, stepVal, windowSize));
    return rcpp_result_gen;
END_RCPP
}
// GetGenomeId
std::string GetGenomeId(std::string inStr);
static SEXP _cognac_GetGenomeId_try(SEXP inStrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< std::string >::type inStr(inStrSEXP);
    rcpp_result_gen = Rcpp::wrap(GetGenomeId(inStr));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _cognac_GetGenomeId(SEXP inStrSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_cognac_GetGenomeId_try(inStrSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// ParseCdHit
void ParseCdHit(const std::string& cdHitClstFile, bool isBinary, int clSizeThesh, Rcpp::Environment& geneEnv);
RcppExport SEXP _cognac_ParseCdHit(SEXP cdHitClstFileSEXP, SEXP isBinarySEXP, SEXP clSizeTheshSEXP, SEXP geneEnvSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type cdHitClstFile(cdHitClstFileSEXP);
    Rcpp::traits::input_parameter< bool >::type isBinary(isBinarySEXP);
    Rcpp::traits::input_parameter< int >::type clSizeThesh(clSizeTheshSEXP);
    Rcpp::traits::input_parameter< Rcpp::Environment& >::type geneEnv(geneEnvSEXP);
    ParseCdHit(cdHitClstFile, isBinary, clSizeThesh, geneEnv);
    return R_NilValue;
END_RCPP
}
// ParseFasta
Rcpp::CharacterVector ParseFasta(const std::string& faPath);
RcppExport SEXP _cognac_ParseFasta(SEXP faPathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type faPath(faPathSEXP);
    rcpp_result_gen = Rcpp::wrap(ParseFasta(faPath));
    return rcpp_result_gen;
END_RCPP
}
// TranslateAaAlgnToDna
void TranslateAaAlgnToDna(const Rcpp::DataFrame& gffData, const std::string& faPath, const std::vector< int >& genePositions, const std::string& genomeName, const std::string& aaAlgn, const std::string& outputFile);
RcppExport SEXP _cognac_TranslateAaAlgnToDna(SEXP gffDataSEXP, SEXP faPathSEXP, SEXP genePositionsSEXP, SEXP genomeNameSEXP, SEXP aaAlgnSEXP, SEXP outputFileSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type gffData(gffDataSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type faPath(faPathSEXP);
    Rcpp::traits::input_parameter< const std::vector< int >& >::type genePositions(genePositionsSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type genomeName(genomeNameSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type aaAlgn(aaAlgnSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type outputFile(outputFileSEXP);
    TranslateAaAlgnToDna(gffData, faPath, genePositions, genomeName, aaAlgn, outputFile);
    return R_NilValue;
END_RCPP
}

// validate (ensure exported C++ functions exist before calling them)
static int _cognac_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("std::string(*GetGenomeId)(std::string)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _cognac_RcppExport_registerCCallable() { 
    R_RegisterCCallable("cognac", "_cognac_GetGenomeId", (DL_FUNC)_cognac_GetGenomeId_try);
    R_RegisterCCallable("cognac", "_cognac_RcppExport_validate", (DL_FUNC)_cognac_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_cognac_CalcAlgnSubMatrix", (DL_FUNC) &_cognac_CalcAlgnSubMatrix, 1},
    {"_cognac_CalcAlgnPartitionDists", (DL_FUNC) &_cognac_CalcAlgnPartitionDists, 3},
    {"_cognac_ConcatenateAlignments", (DL_FUNC) &_cognac_ConcatenateAlignments, 2},
    {"_cognac_CreateAlgnDistMat", (DL_FUNC) &_cognac_CreateAlgnDistMat, 2},
    {"_cognac_CreateCognacRunData", (DL_FUNC) &_cognac_CreateCognacRunData, 4},
    {"_cognac_CreateCoreGenomeDistMat", (DL_FUNC) &_cognac_CreateCoreGenomeDistMat, 1},
    {"_cognac_DeletePartitions", (DL_FUNC) &_cognac_DeletePartitions, 4},
    {"_cognac_ExtractGenomeNameFromPath", (DL_FUNC) &_cognac_ExtractGenomeNameFromPath, 1},
    {"_cognac_GetGenomeNameWithExt", (DL_FUNC) &_cognac_GetGenomeNameWithExt, 2},
    {"_cognac_FilterAlgnPositions", (DL_FUNC) &_cognac_FilterAlgnPositions, 4},
    {"_cognac_FilterPartitionedAlgnPositions", (DL_FUNC) &_cognac_FilterPartitionedAlgnPositions, 5},
    {"_cognac_FindIdenticalGenes", (DL_FUNC) &_cognac_FindIdenticalGenes, 2},
    {"_cognac_GetAlgnQualScores", (DL_FUNC) &_cognac_GetAlgnQualScores, 4},
    {"_cognac_GetGenomeId", (DL_FUNC) &_cognac_GetGenomeId, 1},
    {"_cognac_ParseCdHit", (DL_FUNC) &_cognac_ParseCdHit, 4},
    {"_cognac_ParseFasta", (DL_FUNC) &_cognac_ParseFasta, 1},
    {"_cognac_TranslateAaAlgnToDna", (DL_FUNC) &_cognac_TranslateAaAlgnToDna, 6},
    {"_cognac_RcppExport_registerCCallable", (DL_FUNC) &_cognac_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_cognac(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
