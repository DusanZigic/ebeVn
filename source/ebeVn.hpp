#ifndef HEADERFILE_EBEVNHEADER
#define HEADERFILE_EBEVNHEADER

class ebeVn {
public:
    ebeVn(int argc, const char *argv[]);
    ~ebeVn();
    void run();

private:
    bool m_error = false;

    std::string m_collsys;
    std::string m_sNN;
    std::string m_method;
    std::string m_pName;
    std::string m_centrality;
    std::string m_xB;
    std::vector<int> m_filter;

    size_t m_eventN;
    std::map<size_t, double> m_M;
    std::vector<unsigned int> m_nList;
    std::map<size_t, std::map<unsigned int, std::complex<double>>> m_Qn, m_normedQn;
    std::map<size_t, std::map<unsigned int, double>> m_Psin;
    int loadQnVectors();

    std::vector<double> m_pTGrid, m_phiGrid;
    int setGrids();

    std::map<std::string, std::vector<double>> m_gaussQuadPts;
    int setGaussQuadrature();

    int loadRAADist(int eventID, std::vector<std::vector<double>> &RAADist);

    void integrateRAASP(const std::vector<std::vector<double>> &RAADist, std::vector<double> &RAA, std::map<unsigned int, std::vector<std::complex<double>>> &qn);
    int generateVnSP(std::vector<double> &RAA, std::map<unsigned int, std::vector<double>> &Vn);

    int exportObservables(std::vector<double> &RAA, std::map<unsigned int, std::vector<double>> &Vn);
};

#endif