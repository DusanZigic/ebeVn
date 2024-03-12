#ifndef HEADERFILE_EBEVNHEADER
#define HEADERFILE_EBEVNHEADER

#include <string>
#include <vector>
#include <map>
#include <complex>

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
    double m_xB;
    std::vector<unsigned int> m_filter;

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

    int loadRAADist(size_t eventID, std::vector<std::vector<double>> &RAADist) const;

    void integrateRAASP(const std::vector<std::vector<double>> &RAADist, std::vector<double> &RAA, std::map<unsigned int, std::vector<std::complex<double>>> &qn);
    int calculateVnSP(std::vector<double> &RAA, std::map<unsigned int, std::vector<double>> &Vn);

    void integrateRAAEP(size_t eventID, const std::vector<std::vector<double>> &RAADist, std::vector<double> &RAA, std::map<unsigned int, std::vector<double>> &vn);
    int calculateVnEP(std::vector<double> &RAA, std::map<unsigned int, std::vector<double>> &Vn);

    int exportObservables(const std::vector<double> &RAA, const std::map<unsigned int, std::vector<double>> &Vn) const;

    void integrateRAAC(const std::vector<std::vector<double>> &RAADist, std::vector<double> &mq, std::map<unsigned int, std::vector<std::complex<double>>> &qn);
    void softCummulants(std::map<unsigned int, double> &c2, std::map<unsigned int, double> &C2, std::map<unsigned int, double> &c4, std::map<unsigned int, double> &C4);
    int  hardCummulants(std::vector<double> &RAA, std::map<unsigned int, std::vector<double>> &d2, std::map<unsigned int, std::vector<double>> &d4);
    int  calculateVnC(std::vector<double> &RAA, std::map<unsigned int, std::vector<double>> &Vn2, std::map<unsigned int, std::vector<double>> &Vn4);

    int  exportObservables(const std::vector<double> &RAA, const std::map<unsigned int, std::vector<double>> &Vn2, const std::map<unsigned int, std::vector<double>> &Vn4) const;
};

#endif