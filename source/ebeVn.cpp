#include "ebeVn.hpp"

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <complex>
#include <map>
#include <fstream>
#include <iomanip>

ebeVn::ebeVn(int argc, const char *argv[])
{
    std::vector<std::string> inputs; for (int i=1; i<argc; i++) inputs.push_back(argv[i]);
    
    if ((inputs.size() == 1) && (inputs[0] == "-h")) {
		std::cout << "default values: --collsys=PbPb --sNN=5020GeV --method=SP --filter=off --pName=Charm --centrality=30-40% --xB=0.6" << std::endl;
		m_error = true;
	}

    std::map<std::string, std::string> inputparams;
	for (auto in : inputs)
	{
 	   	std::string key = in.substr(0, in.find("="));
 	   	std::string::size_type n = 0; while ((n = key.find("-", n)) != std::string::npos) {key.replace(n, 1, ""); n += 0;} //replacing all '-'
		std::string val = in.substr(in.find("=")+1, in.length());
		inputparams[key] = val;
	}

    //checking if configuration file is provided:
	std::map<std::string, std::string> inputparamsFile;
	if (inputparams.count("c") > 0) {
		std::ifstream file_in(inputparams["c"], std::ios_base::in);
		if (!file_in.is_open()) {
			std::cerr << "Error: unable to open configuration file. Aborting..." << std::endl;
			m_error = true;
		}
		std::string line, key, sep, val;
		while (std::getline(file_in, line))
		{
			std::stringstream ss(line);
			ss >> key; ss >>sep; ss >> val;
			inputparamsFile[key] = val;
		}
		file_in.close();
	}

    //setting parameter values based on config file values and overwriting with command line values:
	//
	m_collsys = "PbPb"; if (inputparamsFile.count("collsys") > 0) m_collsys = inputparamsFile["collsys"];
						if (    inputparams.count("collsys") > 0) m_collsys =     inputparams["collsys"];
						    
	m_sNN = "5020GeV"; if (inputparamsFile.count("sNN") > 0) m_sNN = inputparamsFile["sNN"];
					   if (    inputparams.count("sNN") > 0) m_sNN =     inputparams["sNN"];

	m_method = "SP"; if (inputparamsFile.count("method") > 0) m_method = inputparamsFile["method"];
					 if (    inputparams.count("method") > 0) m_method =     inputparams["method"];

	std::string filterstr = "0"; if (inputparamsFile.count("filter") > 0) filterstr  =  inputparamsFile["filter"];
							     if (    inputparams.count("filter") > 0) filterstr  =      inputparams["filter"];

	m_pName = "Charm"; if (inputparamsFile.count("pName") > 0) m_pName = inputparamsFile["pName"];
					   if (    inputparams.count("pName") > 0) m_pName =     inputparams["pName"];

	m_centrality = "30-40%"; if (inputparamsFile.count("centrality") > 0) m_centrality = inputparamsFile["centrality"];
						     if (    inputparams.count("centrality") > 0) m_centrality =     inputparams["centrality"];

	m_xB = 0.6; if (inputparamsFile.count("xB") > 0) m_xB = std::stod(inputparamsFile["xB"]);
				if (    inputparams.count("xB") > 0) m_xB = std::stod(    inputparams["xB"]);


	//checking if provided value of sNN is an option:
	if ((m_sNN != "5440GeV") && (m_sNN != "5020GeV") && (m_sNN != "2760GeV") && (m_sNN != "200GeV")) {
		std::cerr << "Error: provided sNN parameter not an option, please try 5440GeV, 5020GeV, 2760GeV or 200GeV. Aborting..." << std::endl;
		m_error = true;
	}

	//setting filter integer values:
	int buffer; std::stringstream ss(filterstr); ss >> buffer;
	if (buffer != 0) {
		while(buffer != 0) {m_filter.push_back(buffer % 10); buffer /= 10;}
	}
	std::reverse(m_filter.begin(), m_filter.end());
}

ebeVn::~ebeVn() {}

void ebeVn::run()
{
    if (loadQnVectors()      != 1) return;
    if (setGrids()           != 1) return;
    if (setGaussQuadrature() != 1) return;

    if (m_method == "SP") {
        std::vector<double> RAA; std::map<unsigned int, std::vector<double>> Vn;
        if (     generateVnSP(RAA, Vn) != 1) return;
        if (exportObservables(RAA, Vn) != 1) return;
    }
}

int ebeVn::loadQnVectors()
{
    std::string path_in = "../QnVectors/QnVectors_cent=" + m_centrality + ".dat";
	std::ifstream file_in(path_in, std::ios_base::in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open Qn vectors file. Aborting..." << std::endl;
		return -1;
	}

	std::string line, bufferS; double buffM, buffRe, buffIm; size_t eventID;

	while (getline(file_in, line)) {
		if (line.at(0) == '#') {
            std::stringstream ss(line);
            while (ss >> bufferS) {
                if (bufferS.find("Q") != std::string::npos)
                    m_nList.push_back(std::stoi(bufferS.substr(0, bufferS.find("_")).substr(1, 1)));
            }
            std::sort(m_nList.begin(), m_nList.end());
            m_nList.erase(std::unique(m_nList.begin(), m_nList.end()), m_nList.end());
			continue;
        }

		std::stringstream ss(line);

		ss >> eventID;

		ss >> buffM; m_M[eventID] = buffM;

        for (const auto& n : m_nList) {
            ss >> buffRe;
            ss >> buffIm;
            m_Qn[eventID][n] = std::complex<double>{buffRe, buffIm};
        }

	}

	file_in.close();

    m_eventN = m_M.size();

    for (const auto &qn : m_Qn) {
        eventID = qn.first;
        for (const auto& n : m_nList) {
            m_normedQn[eventID][n] = m_Qn[eventID][n] / m_M[eventID];
                m_Psin[eventID][n] = 1.0/static_cast<double>(n)*std::atan2(m_normedQn[eventID][n].imag(), m_normedQn[eventID][n].real());
        }
    }
}

int ebeVn::setGrids()
{
    std::stringstream xBsstr; xBsstr << std::fixed << std::setprecision(1) << m_xB;
    size_t eventID = m_Qn.begin()->first;
    std::string path_in  = "../results/results" + m_pName + "/";
                path_in += m_pName + "_" + m_collsys + "_sNN=" + m_sNN + "_cent=" + m_centrality + "_xB=" + xBsstr.str();
                path_in += "_dist_" + std::to_string(eventID) + ".dat";
    std::ifstream file_in(path_in, std::ios_base::in);
    if (!file_in.is_open()) {
        std::cerr << "Error: unable to open R_AA(pT,phi) distribution file while generating grids. Aborting..." << std::endl;
        return -1;
    }

    std::string line; double buffer;

    while (std::getline(file_in, line)) {
        if (line.at(0) == '#')
            continue;
        std::stringstream ss(line);
        ss >> buffer; m_pTGrid.push_back(buffer);
        ss >> buffer; m_phiGrid.push_back(buffer);
        ss >> buffer;
    }

    file_in.close();

    std::sort(m_pTGrid.begin(),  m_pTGrid.end());  m_pTGrid.erase(std::unique( m_pTGrid.begin(),  m_pTGrid.end()),  m_pTGrid.end());
    std::sort(m_phiGrid.begin(), m_phiGrid.end()); m_phiGrid.erase(std::unique(m_phiGrid.begin(), m_phiGrid.end()), m_phiGrid.end());

    return 1;
}

int ebeVn::setGaussQuadrature()
{
    std::string path_in = "./phiGaussPts/phiptsgauss" + std::to_string(m_phiGrid.size()) + ".dat";
    std::ifstream file_in(path_in, std::ios_base::in);
    if (!file_in.is_open()) {
        std::cerr << "Error: unable to open Gaussian quadtrature points and weights file. Aborting..." << std::endl;
        return -1;
    }

    m_gaussQuadPts["phi"]     = std::vector<double>();
    m_gaussQuadPts["weights"] = std::vector<double>();

    std::string line; double buffer;

    while (std::getline(file_in, line)) {
        if (line.at(0) == '#')
            continue;
        std::stringstream ss(line);
        ss >> buffer; m_gaussQuadPts["phi"].push_back(buffer);
        ss >> buffer; m_gaussQuadPts["weights"].push_back(buffer);
    }    

    file_in.close();

    if (m_gaussQuadPts["phi"].size() != m_phiGrid.size()) {
        std::cerr << "Error: Gaussian quadrature phi points and phi grids are not the same size. Aborting..." << std::endl;
        return -2;
    }

    return 1;
}

int ebeVn::loadRAADist(int eventID, std::vector<std::vector<double>> &RAADist)
{
	std::stringstream xBsstr; xBsstr << std::fixed << std::setprecision(1) << m_xB;
    std::string path_in  = "../results/results" + m_pName + "/";
                path_in += m_pName + "_" + m_collsys + "_sNN=" + m_sNN + "_cent=" + m_centrality + "_xB=" + xBsstr.str();
                path_in += "_dist_" + std::to_string(eventID) + ".dat";
    std::ifstream file_in(path_in, std::ios_base::in);
    if (!file_in.is_open()) {
        std::cerr << "Error: unable to open R_AA(pT,phi) distribution file for event: " + std::to_string(eventID) + ". Aborting..." << std::endl;
        return -1;
    }

	std::vector<double> pTPts, phiPts, RAAPts;
	std::string line; double buffer;

	while(std::getline(file_in, line))
	{
		if (line.at(0) == '#')
            continue;

		std::stringstream ss(line);
		ss >> buffer;  pTPts.push_back(buffer);
		ss >> buffer; phiPts.push_back(buffer);
		ss >> buffer; RAAPts.push_back(buffer);
	}

    RAADist.resize(m_pTGrid.size(), std::vector<double>(m_phiGrid.size(), 0.0));

	for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++)
		for (size_t iphi=0; iphi<m_phiGrid.size(); iphi++)
			RAADist[ipT][iphi] = RAAPts[ipT*m_phiGrid.size() + iphi];

	file_in.close();

	return 1;
}


void ebeVn::integrateRAASP(const std::vector<std::vector<double>> &RAADist, std::vector<double> &RAA, std::map<unsigned int, std::vector<std::complex<double>>> &qn)
{
    std::vector<double> norm(m_pTGrid.size(), 0.0);
	for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++)
		for (size_t iphi=0; iphi<m_gaussQuadPts["phi"].size(); iphi++)
			norm[ipT] += m_gaussQuadPts["weights"][iphi]*RAADist[ipT][iphi];

	RAA.resize(m_pTGrid.size(), 0.0);
	for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++)
		RAA[ipT] = norm[ipT]/2.0/M_PI;
    
    for (const auto &n : m_nList) {
        qn[n] = std::vector<std::complex<double>>(m_pTGrid.size(), 0.0);
        for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++) {
		    for (size_t iphi=0; iphi<m_gaussQuadPts["phi"].size(); iphi++) {
                qn[n][ipT] += m_gaussQuadPts["weights"][iphi]*RAADist[ipT][iphi]*std::exp(std::complex<double>{0.0, static_cast<double>(n)*m_phiGrid[iphi]});
            }
            qn[n][ipT] /= norm[ipT];
        }
    }
}

int ebeVn::generateVnSP(std::vector<double> &RAA, std::map<unsigned int, std::vector<double>> &Vn)
{
    return 1;
}

int ebeVn::exportObservables(std::vector<double> &RAA, std::map<unsigned int, std::vector<double>> &Vn)
{
    return 1;
}