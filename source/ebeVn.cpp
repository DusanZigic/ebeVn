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
	std::sort(m_filter.begin(), m_filter.end());
}

ebeVn::~ebeVn() {}

void ebeVn::run()
{
    if (loadQnVectors()      != 1) return;
    if (setGrids()           != 1) return;
    if (setGaussQuadrature() != 1) return;

    if (m_method == "SP") {
        std::vector<double> RAA; std::map<unsigned int, std::vector<double>> Vn;
        if (    calculateVnSP(RAA, Vn) != 1) return;
        if (exportObservables(RAA, Vn) != 1) return;
    }
    else if (m_method == "EP") {
        std::vector<double> RAA; std::map<unsigned int, std::vector<double>> Vn;
        if (    calculateVnEP(RAA, Vn) != 1) return;
        if (exportObservables(RAA, Vn) != 1) return;
    }
    else if (m_method == "C") {
        std::vector<double> RAA; std::map<unsigned int, std::vector<double>> Vn2, Vn4;
        if (     calculateVnC(RAA, Vn2, Vn4) != 1) return;
        if (exportObservables(RAA, Vn2, Vn4) != 1) return;
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
            for (const auto &n : m_nList) {
				if (n == 1) continue;
				if(std::find(m_nList.begin(), m_nList.end(), 2*n) == m_nList.end()) continue;
				m_nListCummulants.push_back(n);
			}
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

    return 1;
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
        ss >> buffer; m_gaussQuadPts.at("phi").push_back(buffer);
        ss >> buffer; m_gaussQuadPts.at("weights").push_back(buffer);
    }    

    file_in.close();

    if (m_gaussQuadPts.at("phi").size() != m_phiGrid.size()) {
        std::cerr << "Error: Gaussian quadrature phi points and phi grids are not the same size. Aborting..." << std::endl;
        return -2;
    }

    return 1;
}

int ebeVn::loadRAADist(size_t eventID, std::vector<std::vector<double>> &RAADist) const
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
		for (size_t iphi=0; iphi<m_gaussQuadPts.at("phi").size(); iphi++)
			norm[ipT] += m_gaussQuadPts.at("weights")[iphi]*RAADist[ipT][iphi];

	RAA.resize(m_pTGrid.size(), 0.0);
	for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++)
		RAA[ipT] = norm[ipT]/2.0/M_PI;
    
    for (const auto &n : m_nList) {
        qn[n] = std::vector<std::complex<double>>(m_pTGrid.size(), std::complex<double>{0.0, 0.0});
        for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++) {
		    for (size_t iphi=0; iphi<m_gaussQuadPts.at("phi").size(); iphi++) {
                qn[n][ipT] += m_gaussQuadPts.at("weights")[iphi]*RAADist[ipT][iphi]*std::exp(std::complex<double>{0.0, static_cast<double>(n)*m_phiGrid[iphi]});
            }
            qn[n][ipT] /= norm[ipT];
        }
    }
}

int ebeVn::calculateVnSP(std::vector<double> &RAA, std::map<unsigned int, std::vector<double>> &Vn)
{
    RAA.resize(m_pTGrid.size(), 0.0);
    for (const auto &n : m_nList)
        Vn[n] = std::vector<double>(m_pTGrid.size(), 0.0);

    std::map<unsigned int, std::vector<double>> Qqn;
    std::map<unsigned int, double> QQn;
    for (const auto &n : m_nList) {
        Qqn[n] = std::vector<double>(m_pTGrid.size(), 0.0);
        QQn[n] = 0.0;
    }    
 
    size_t filterEventN = 0;

    for (const auto & Qn : m_normedQn) {
        size_t eventID = Qn.first;
        
        bool skipEvent = false;
		for (const auto &f: m_filter) {
            if (std::abs(m_Qn[eventID][2*f]) > std::abs(m_Qn[eventID][f])) {
                skipEvent = true;
            }
        }
		if (skipEvent)
            continue;
        
        filterEventN++;

        std::vector<std::vector<double>> RAADist;
        if (loadRAADist(eventID, RAADist) != 1) return -1;

		std::vector<double> raa; std::map<unsigned int, std::vector<std::complex<double>>> qn;
        integrateRAASP(RAADist, raa, qn);

		for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++) {
			RAA[ipT] += raa[ipT];
			for (const auto &n : m_nList)
				Qqn[n][ipT] += std::real(qn[n][ipT]*std::conj(m_normedQn[eventID][n]));
		}

		for (const auto &n : m_nList)
			QQn[n] += std::real(m_normedQn[eventID][n]*std::conj(m_normedQn[eventID][n]));
    }

    for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++)
	{
		RAA[ipT] /= static_cast<double>(filterEventN);
		for (const auto &n : m_nList)
			Vn[n][ipT] = std::real(Qqn[n][ipT])/static_cast<double>(filterEventN) / std::sqrt(QQn[n]/static_cast<double>(filterEventN));
	}

    return 1;
}


void ebeVn::integrateRAAEP(size_t eventID, const std::vector<std::vector<double>> &RAADist, std::vector<double> &RAA, std::map<unsigned int, std::vector<double>> &vn)
{
    std::vector<double> norm(m_pTGrid.size(), 0.0);
	for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++)
		for (size_t iphi=0; iphi<m_phiGrid.size(); iphi++)
			norm[ipT] += m_gaussQuadPts.at("weights")[iphi]*RAADist[ipT][iphi];

    RAA.resize(m_pTGrid.size(), 0.0);
	for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++)
		RAA[ipT] = norm[ipT]/2.0/M_PI;

	for (const auto &n : m_nList) {
		vn[n] = std::vector<double>(m_pTGrid.size(), 0.0);

		for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++) {
			for (size_t iphi=0; iphi<m_gaussQuadPts.at("phi").size(); iphi++)
				vn[n][ipT] += m_gaussQuadPts.at("weights")[iphi]*RAADist[ipT][iphi]*std::cos(static_cast<double>(n)*(m_gaussQuadPts.at("phi")[iphi] - m_Psin[eventID][n]));

			vn[n][ipT] /= norm[ipT];
		}
	}
}

int ebeVn::calculateVnEP(std::vector<double> &RAA, std::map<unsigned int, std::vector<double>> &Vn)
{
    RAA.resize(m_pTGrid.size(), 0.0);
    for (const auto &n : m_nList)
        Vn[n] = std::vector<double>(m_pTGrid.size(), 0.0);
    
    size_t filterEventN = 0;

    for (const auto & Psin : m_Psin) {
        size_t eventID = Psin.first;
        
        bool skipEvent = false;
		for (const auto &f: m_filter) {
            if (std::abs(m_Qn[eventID][2*f]) > std::abs(m_Qn[eventID][f])) {
                skipEvent = true;
            }
        }
		if (skipEvent)
            continue;
        
        filterEventN++;

        std::vector<std::vector<double>> RAADist;
        if (loadRAADist(eventID, RAADist) != 1) return -1;

        std::vector<double> raa; std::map<unsigned int, std::vector<double>> vn;
        integrateRAAEP(eventID, RAADist, raa, vn);

        for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++) {
            RAA[ipT] += raa[ipT];
            for (const auto &n : m_nList)
                Vn[n][ipT] += vn[n][ipT];
        }
    }

    for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++) {
        RAA[ipT] /= static_cast<double>(filterEventN);
        for (const auto &n : m_nList)
            Vn[n][ipT] /= static_cast<double>(filterEventN);
    }

    return 1;
}

int ebeVn::exportObservables(const std::vector<double> &RAA, const std::map<unsigned int, std::vector<double>> &Vn) const
{
    std::stringstream xBsstr; xBsstr << std::fixed << std::setprecision(1) << m_xB;

	std::string path_out  = "../results/results" + m_pName + "/";
	            path_out += m_pName + "_" + m_collsys + "_sNN=" + m_sNN + "_cent=" + m_centrality + "_xB=" + xBsstr.str();

	if      (m_method == "SP") path_out += "_SP";
	else if (m_method == "EP") path_out += "_EP";

	if (m_filter.size() > 0) {
		path_out += "_F";
        for (const auto &f : m_filter)
            path_out = path_out + std::to_string(f);
	}

	path_out += ".dat";

	std::ofstream file_out(path_out, std::ios_base::out);
	if (!file_out.is_open()) {
		std::cerr << "Error: unable top open output file. Aborting..." << std::endl;
		return -1;
	}

    file_out << "#";
    file_out << std::fixed << std::setw(13) << "pT [GeV]" << " ";
    file_out << std::fixed << std::setw(13) <<     "R_AA" << " ";
	if (m_method == "SP") {
		for (const auto &n : m_nList) {
            std::string vnStr = "v_" + std::to_string(n) + "{SP}";
            file_out << std::fixed << std::setw(13) << vnStr << " ";
        }
	}
	else if (m_method == "EP") {
		for (const auto &n : m_nList) {
            std::string vnStr = "v_" + std::to_string(n) + "{EP}";
            file_out << std::fixed << std::setw(13) << vnStr << " ";
        }
	}
    file_out << "\n";

	for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++) {
		file_out << std::fixed << std::setw(14) << std::setprecision(10) << m_pTGrid[ipT] << " ";
		file_out << std::fixed << std::setw(13) << std::setprecision(10) << 	 RAA[ipT] << " ";
		for (const auto &n : m_nList)
		    file_out << std::fixed << std::setw(13) << std::setprecision(10) <<  Vn.at(n)[ipT] << " ";
		file_out << "\n";
	}

	return 1;
}


void ebeVn::integrateRAAC(const std::vector<std::vector<double>> &RAADist, std::vector<double> &mq, std::map<unsigned int, std::vector<std::complex<double>>> &qn)
{
    mq.resize(m_pTGrid.size(), 0.0);
	for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++)
		for (size_t iphi=0; iphi<m_gaussQuadPts.at("phi").size(); iphi++)
			mq[ipT] += m_gaussQuadPts.at("weights")[iphi]*RAADist[ipT][iphi];

	for (const auto &n : m_nList) {
		qn[n] = std::vector<std::complex<double>>(m_pTGrid.size(), std::complex<double>{0.0, 0.0});

		for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++) {
			for (size_t iphi=0; iphi<m_gaussQuadPts.at("phi").size(); iphi++)
				qn[n][ipT] += m_gaussQuadPts.at("weights")[iphi]*RAADist[ipT][iphi]*std::exp(std::complex<double>{0.0, static_cast<double>(n)*m_phiGrid[iphi]});
		}
	}
}

void ebeVn::softCummulants(std::map<unsigned int, double> &c2, std::map<unsigned int, double> &C2, std::map<unsigned int, double> &c4, std::map<unsigned int, double> &C4)
{
    //legend: c2 = <<2>>, C2 = c_n{2}, c4 = <<4>>, C4 = c_n{4}

	std::map<unsigned int, double>  twoParticleCummulantSum; double  twoParticleCummulantNorm = 0.0;
	std::map<unsigned int, double> fourParticleCummulantSum; double fourParticleCummulantNorm = 0.0;

    for (const auto& n : m_nListCummulants) {
         twoParticleCummulantSum[n] = 0.0;
        fourParticleCummulantSum[n] = 0.0;
    }

    size_t filterEventN = 0;

	for (const auto & Qn : m_Qn) {
        size_t eventID = Qn.first;
        
        bool skipEvent = false;
		for (const auto &f: m_filter) {
            if (std::abs(m_Qn[eventID][2*f]) > std::abs(m_Qn[eventID][f])) {
                skipEvent = true;
            }
        }
		if (skipEvent)
            continue;
        
        filterEventN++;

		double W2 = m_M[eventID]*(m_M[eventID]-1.0);
		double W4 = m_M[eventID]*(m_M[eventID]-1.0)*(m_M[eventID]-2.0)*(m_M[eventID]-3.0);

		for (const auto& n : m_nListCummulants) {
		     twoParticleCummulantSum[n] += ((std::pow(std::abs(m_Qn[eventID][n]), 2.0) - m_M[eventID])/W2) * W2;
		    fourParticleCummulantSum[n] += ((std::pow(std::abs(m_Qn[eventID][n]), 4.0) + std::pow(std::abs(m_Qn[eventID][2*n]), 2.0) - 
                                                2.0*std::real(m_Qn[eventID][2*n]*std::conj(m_Qn[eventID][n])*std::conj(m_Qn[eventID][n])))/W4 -
                                            2.0*(2.0*(m_M[eventID]-2.0)*std::pow(std::abs(m_Qn[eventID][n]), 2.0) -
                                                m_M[eventID]*(m_M[eventID]-3.0))/W4) * W4;
		}

		 twoParticleCummulantNorm += W2;
        fourParticleCummulantNorm += W4;
	}

	for (const auto& n : m_nListCummulants) {
		c2[n] =  twoParticleCummulantSum[n]/ twoParticleCummulantNorm;
        c4[n] = fourParticleCummulantSum[n]/fourParticleCummulantNorm;
	}

	for (const auto& n : m_nListCummulants) {
		C2[n] = c2[n];
		C4[n] = c4[n] - 2.0*c2[n]*c2[n];
	}
}

int ebeVn::hardCummulants(std::vector<double> &RAA, std::map<unsigned int, std::vector<double>> &d2, std::map<unsigned int, std::vector<double>> &d4)
{
	//legend: d2 = <2'>, D2 = d_n{2}, d4 = <4'>, D4 = d_n{4}

    std::map<unsigned int, std::vector<double>>  twoParticleCummulantSum; std::vector<double>  twoParticleCummulantNorm(m_pTGrid.size(), 0.0);
	std::map<unsigned int, std::vector<double>> fourParticleCummulantSum; std::vector<double> fourParticleCummulantNorm(m_pTGrid.size(), 0.0);

    for (const auto& n : m_nListCummulants) {
         twoParticleCummulantSum[n] = std::vector<double>(m_pTGrid.size(), 0.0);
        fourParticleCummulantSum[n] = std::vector<double>(m_pTGrid.size(), 0.0);
    }

	size_t filterEventN = 0;

    RAA.resize(m_pTGrid.size(), 0.0);

	for (const auto & Qn : m_Qn) {
		size_t eventID = Qn.first;
        
        bool skipEvent = false;
		for (const auto &f: m_filter) {
            if (std::abs(m_Qn[eventID][2*f]) > std::abs(m_Qn[eventID][f])) {
                skipEvent = true;
            }
        }
		if (skipEvent)
            continue;
        
        filterEventN++;

        std::vector<std::vector<double>> RAADist;
        if (loadRAADist(eventID, RAADist) != 1) return -1;

        std::vector<double> mq; std::map<unsigned int, std::vector<std::complex<double>>> qn;
        integrateRAAC(RAADist, mq, qn);


		for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++)
            RAA[ipT] += (mq[ipT]/2.0/M_PI);

		std::vector<double> W2(m_pTGrid.size(), 0.0), W4(m_pTGrid.size(), 0.0);

		for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++) {
			W2[ipT] = mq[ipT]*m_M[eventID];
			W4[ipT] = mq[ipT]*m_M[eventID]*(m_M[eventID]-1.0)*(m_M[eventID]-2.0);
			
			for (const auto& n : m_nListCummulants) {
				 twoParticleCummulantSum[n][ipT] += std::real(qn[n][ipT]*std::conj(m_Qn[eventID][n])/W2[ipT]) * W2[ipT];
				fourParticleCummulantSum[n][ipT] += std::real(qn[n][ipT]*(m_Qn[eventID][n]*std::conj(m_Qn[eventID][n])*std::conj(m_Qn[eventID][n]) -
                                                                            m_Qn[eventID][n]*std::conj(m_Qn[eventID][2*n]) -
														                    2.0*m_M[eventID]*std::conj(m_Qn[eventID][n]) + 2.0*std::conj(m_Qn[eventID][n]))/W4[ipT])*W4[ipT];
			}

			 twoParticleCummulantNorm[ipT] += W2[ipT];
            fourParticleCummulantNorm[ipT] += W4[ipT];
		}
	}

	for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++)
        RAA[ipT] /= static_cast<double>(filterEventN);

	for (const auto& n : m_nListCummulants) {
        d2[n] = std::vector<double>(m_pTGrid.size(), 0.0);
        d4[n] = std::vector<double>(m_pTGrid.size(), 0.0);
		for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++) {
			d2[n][ipT] =  twoParticleCummulantSum[n][ipT]/ twoParticleCummulantNorm[ipT];
			d4[n][ipT] = fourParticleCummulantSum[n][ipT]/fourParticleCummulantNorm[ipT];
		}
    }

    return 1;
}

int ebeVn::calculateVnC(std::vector<double> &RAA, std::map<unsigned int, std::vector<double>> &Vn2, std::map<unsigned int, std::vector<double>> &Vn4)
{
    if (m_nListCummulants.empty()) {
        std::cerr << "Error: there are no appropriate Qn vectors to calculate cummulants. Aborting..." << std::endl;
        return -1;
    }
    // legend: c2 = <<2>>, C2 = c_n{2}, c4 = <<4>>, C4 = c_n{4}
	std::map<unsigned int, double> c2, C2, c4, C4;
	softCummulants(c2, C2, c4, C4);

    // legend: d2 = <2'>,   d4 = <4'>
	std::map<unsigned int, std::vector<double>> d2, d4;
	if (hardCummulants(RAA, d2, d4) != 1) return -2;

    // legend: //D2 = d_n{2}, D4 = d_n{4}
	std::map<unsigned int, std::vector<double>> D2, D4;
	for (const auto& n : m_nListCummulants) {
        D2[n] = std::vector<double>(m_pTGrid.size(), 0.0);
        D4[n] = std::vector<double>(m_pTGrid.size(), 0.0);
		for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++) {
			D2[n][ipT] = d2[n][ipT];
			D4[n][ipT] = d4[n][ipT] - 2.0*d2[n][ipT]*c2[n];
		}
	}

	for (const auto& n : m_nListCummulants) {
        Vn2[n] = std::vector<double>(m_pTGrid.size(), 0.0);
        Vn4[n] = std::vector<double>(m_pTGrid.size(), 0.0);
		for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++) {
			Vn2[n][ipT] =      D2[n][ipT]/std::pow(     C2[n], 1.0/2.0);
			Vn4[n][ipT] = -1.0*D4[n][ipT]/std::pow(-1.0*C4[n], 3.0/4.0);
		}
	}
    
	return 1;
}

int ebeVn::exportObservables(const std::vector<double> &RAA, const std::map<unsigned int, std::vector<double>> &Vn2, const std::map<unsigned int, std::vector<double>> &Vn4) const
{
    std::stringstream xBsstr; xBsstr << std::fixed << std::setprecision(1) << m_xB;

	std::string path_out  = "../results/results" + m_pName + "/";
	            path_out += m_pName + "_" + m_collsys + "_sNN=" + m_sNN + "_cent=" + m_centrality + "_xB=" + xBsstr.str() + "_C";
	
    if (m_filter.size() > 0) {
		path_out += "_F";
        for (const auto &f : m_filter)
            path_out = path_out + std::to_string(f);
	}

	path_out += ".dat";

	std::ofstream file_out(path_out, std::ios_base::out);
	if (!file_out.is_open()) {
		std::cerr << "Error: unable top open output file. Aborting..." << std::endl;
		return -1;
	}

    file_out << "#";
    file_out << std::fixed << std::setw(13) << "pT [GeV]" << " ";
    file_out << std::fixed << std::setw(13) <<     "R_AA" << " ";
    for (const auto& n : m_nListCummulants) {
        std::string vnStr;
        vnStr = "v_" + std::to_string(n) + "{2}";
        file_out << std::fixed << std::setw(13) << vnStr << " ";
        vnStr = "v_" + std::to_string(n) + "{4}";
        file_out << std::fixed << std::setw(13) << vnStr << " ";
    }
    file_out << "\n";

    for (size_t ipT=0; ipT<m_pTGrid.size(); ipT++) {
		file_out << std::fixed << std::setw(14) << std::setprecision(10) << m_pTGrid[ipT] << " ";
		file_out << std::fixed << std::setw(13) << std::setprecision(10) << 	 RAA[ipT] << " ";
		for (const auto& n : m_nListCummulants) {
		    file_out << std::fixed << std::setw(13) << std::setprecision(10) <<  Vn2.at(n)[ipT] << " ";
            file_out << std::fixed << std::setw(13) << std::setprecision(10) <<  Vn4.at(n)[ipT] << " ";
        }
		file_out << "\n";
	}

	file_out.close();

	return 1;
}