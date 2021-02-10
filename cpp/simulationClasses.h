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

#ifndef _SIMULATIONCLASSES_H
#define _SIMULATIONCLASSES_H
 
#include <vector>
#include <eigen-3.3.9/Eigen/Dense>

using namespace std;
using namespace Eigen;

class sim_output {
public:
    // Constructor:
    sim_output(int l, int n, int m) : length(l), Nx(n), Nu(m) {}; //}, X(n,l), err(n,l), Unom(m,l), Ucorr(m,l), optVal(l), Supervision(l), activeSet(l)  {};
    // Assign values:
    void add_X(VectorXd state) { X.push_back(state); };
    void add_err(VectorXd error) { err.push_back(error); };
    void add_Unom(VectorXd input) { Unom.push_back(input); };
    void add_Ucorr(VectorXd input) { Ucorr.push_back(input); };
    void add_optVal(double val) { optVal.push_back(val); };
    void add_Supervision(bool val) { Supervision.push_back(val); };
    void add_activeSet(int val) { activeSet.push_back(val); };
    // Get values:
    int get_len() { return length; };
    int get_Nx() { return Nx; };
    int get_Nu() { return Nu; };
    vector<VectorXd> get_X(){ return X; };
    vector<VectorXd> get_err(){ return err; };
    vector<VectorXd> get_Unom(){ return Unom; };
    vector<VectorXd> get_Ucorr(){ return Ucorr; };
    vector<double> get_optVal(){ return optVal; };
    vector<bool> get_Supervision(){ return Supervision; };
    vector<int> get_activeSet(){ return activeSet; };

private:
    int length;
    int Nx;
    int Nu;
    vector<VectorXd> X;
    vector<VectorXd> err;
    vector<VectorXd> Unom;
    vector<VectorXd> Ucorr;
    vector<double> optVal;
    vector<bool> Supervision;
    vector<int> activeSet;

};


class trajectory {
public:
    // Constructor:
    trajectory(int n, int m, double dt, int l, VectorXd gain) : Nx(n), Nu(m), dt(dt), length(l), K(gain), state(Nx,length), input(Nu,length) {};
    // Assign values:
    void add_state(VectorXd new_state, int idx) { state.col(idx) = new_state; };
    void add_input(VectorXd new_input, int idx) { input.col(idx) = new_input; };
    // Get values:
    VectorXd get_state(int idx){ return state.col(idx); };
    VectorXd get_input(int idx){ return input.col(idx); };
    // Get u_des:
    VectorXd nominalInput(int idx, VectorXd error){
        int m = input.rows();
        int n = error.rows()/3;
        VectorXd u_des = VectorXd(m);
        for (int i = 0; i < m; i++)
            u_des[i] = K(0)*error(i) + K(1)*error(i+n) + K(2)*error(i+2*n) + input(i,idx);
        return u_des;
    }
private:
    int Nx, Nu;
    double dt;
    int length;
    VectorXd K;
    MatrixXd state;
    MatrixXd input;
};

struct model {
    int Nx;
    int Nu;
    MatrixXd Ad;
    MatrixXd Bd;
    // Constructor:
    model(int n, int m, MatrixXd A, MatrixXd B) : Nx(n), Nu(m), Ad(A), Bd(B) {};
};

struct polytope {
    MatrixXd A;
    VectorXd b;
    // Constructor:
    polytope(MatrixXd A, MatrixXd b) : A(A), b(b) {};
};

#endif
