/*  Copyright (C) 2020, Tzanis Anevlavis.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.  */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <fstream>
#include "gurobiRoutines.h"
#include "simulationClasses.h"
#include "supervisionRoutines.h"
#include </Users/j10/prog-libs/eigen-3.3.9/Eigen/Dense>

using namespace std;
using namespace Eigen;

// Write Eigen::MatrixXd to CSV:
const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");
void writeToCSVfile(string name, MatrixXd matrix){
    ofstream myfile(name.c_str());
    myfile << matrix.format(CSVFormat);
    myfile.close();
}
// Read matrix from CSV to Eigen:
template<typename M>
M load_csv (const string & path) {
    ifstream indata;
    indata.open(path);
    string line;
    vector<double> values;
    uint rows = 0;
    while (getline(indata, line)) {
        stringstream lineStream(line);
        string cell;
        while (getline(lineStream, cell, ',')) {
            values.push_back(stod(cell));
        }
        ++rows;
    }
    return Map<const Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, RowMajor>>(values.data(), rows, values.size()/rows);
}

int main(){
    // ///////////////////////////////////////////////////////////////////////////////////////////////////
    // LOAD DATA FROM FILES:
    // Load model:
    MatrixXd Ad, Bd;
    Ad = load_csv<MatrixXd>("/Users/j10/Downloads/CIS/quadrocopter_supervision/cpp/gurobi_test/data/Ad.csv");
    Bd = load_csv<MatrixXd>("/Users/j10/Downloads/CIS/quadrocopter_supervision/cpp/gurobi_test/data/Bd.csv");
    model mdl(Ad.cols(), Bd.cols(), Ad, Bd);
    // Create trajectory struct:
    MatrixXd traj_state, traj_input;
    double dt = 0.05;
    int NSim = 1000;
    VectorXd K(3);
    K << 58.191666446869746, 43.645601896558645, 10.883765107576760;
    // Load trajectory:
    traj_state = load_csv<MatrixXd>("/Users/j10/Downloads/CIS/quadrocopter_supervision/cpp/gurobi_test/data/traj_state.csv");
    traj_input = load_csv<MatrixXd>("/Users/j10/Downloads/CIS/quadrocopter_supervision/cpp/gurobi_test/data/traj_input.csv");
    trajectory traj(mdl.Nx, mdl.Nu, dt, NSim, K);
    // Populate trajectory struct:
    for (int i=0; i<NSim; i++){
        traj.add_state(traj_state.col(i),i);
        traj.add_input(traj_input.col(i),i);
    }
    // Load CISs:
    vector<polytope> CISs;
    MatrixXd cisA, cisb;
    for (int i=0; i<9; i++){
        string filename = "/Users/j10/Downloads/CIS/quadrocopter_supervision/cpp/gurobi_test/data/cis"+to_string(i);
        cisA = load_csv<MatrixXd>(filename+"_A.csv");
        cisb = load_csv<MatrixXd>(filename+"_b.csv");
        polytope CIS(cisA, cisb);
        CISs.push_back(CIS);
    }
    // Load input constraints:
    MatrixXd Gu, Fu;
    Gu = load_csv<MatrixXd>("/Users/j10/Downloads/CIS/quadrocopter_supervision/cpp/gurobi_test/data/inputA.csv");
    Fu = load_csv<MatrixXd>("/Users/j10/Downloads/CIS/quadrocopter_supervision/cpp/gurobi_test/data/inputb.csv");
    polytope inputConstr(Gu,Fu);
    // Initial state:
    VectorXd x0(9);
    x0 << -2.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    // ///////////////////////////////////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////////////////////////////////
    // SIMULATION:
    
    // Create simulation struct:
    sim_output sim(NSim, mdl.Nx, mdl.Nu);
    // Initialize state:
    VectorXd x_curr = x0;
    // Start simulation:
    for (int step = 0; step < sim.get_len(); step++){
        // Get trajectory:
        VectorXd x_ref = traj.get_state(step);              // reference next state.
        sim.add_X(x_curr);                                  // store true current state.
        VectorXd err = x_ref - x_curr;                      // trajectory error.
        sim.add_err(err);                                   // store trajectory error.
        VectorXd u_des = traj.nominalInput(step, err);      // nominal input.
        sim.add_Unom(u_des);                                // store nominal input.
        
        VectorXd x_next_des = mdl.Ad*x_curr + mdl.Bd*u_des; //  nominal next state
        // Check if the nominal next state belongs to the CIS.
        // If desired next state belongs to the CIS we are good:
        if (isContained(x_next_des,CISs)){
            x_curr = x_next_des;
            sim.add_Ucorr(u_des);
        }
        // Else we correct the nominal input:
        else{
            cout << "Correcting.. step: " << step << "." << endl;
            sim.add_Supervision(step);
            // We solve for the union of CISs by checking if there exists a next state in any of the CISs separately:
            int Nsets = CISs.size();
            MatrixXd u_cand(mdl.Nu, Nsets); // candidate input
            VectorXd f_cand(Nsets);         // candidate optimal cost
            // For each CIS..
            for (int k = 0; k < Nsets; k++){
                //..call the supervisor:
                pair<double, VectorXd> res = callSupervisor(x_curr, u_des, mdl, CISs[k], inputConstr);
                f_cand(k) = res.first;
                u_cand.col(k) = res.second;
            }

            // Sanity check -- we must have at least one solution:
            if (f_cand.minCoeff()==__DBL_MAX__){
                cout << "Step: " << step << ", out of: " << sim.get_len() << endl;
                cout << "All the optimizations failed.." << endl;
                exit(-1);
            }

            // Pick the input that gives the minimum cost:
            VectorXd::Index next_cis;
            double fval = f_cand.minCoeff(&next_cis);
            sim.add_optVal(fval);
            sim.add_activeSet(next_cis);
            // Select input corresponding to minimum cost:
            VectorXd u_corrected = u_cand.col(next_cis);
            sim.add_Ucorr(u_corrected);
            
            // Make sure that next state is indeed in the next cis:
            VectorXd x_next = mdl.Ad * x_curr + mdl.Bd * u_corrected;
            if (!isContained(x_next, CISs)){
                cout << "Step: " << step << ", out of: " << sim.get_len() << endl;
                //            cisA = DSets{next_cis}.CIS.A;
                //            cisb = DSets{next_cis}.CIS.b;
                //            sat = (cisA * x_next <= cisb);
                //            idcs = find(sat==0);
                //            diff = abs(cisA(idcs,:) * x_next - cisb(idcs));
                cout << "Next state not in a CIS!" << endl;
                exit(-1);
            }
            // Update state:
            x_curr = x_next;
        }

    }
    
    vector<VectorXd> X = sim.get_X();
    cout << X[0] << endl;
    cout << X.back() << endl;
//    writeToCSVfile("fulltrajectory.csv", X);

    
    return 0;
}
