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

#ifndef _SUPERVISIONROUTINES_H
#define _SUPERVISIONROUTINES_H
 
#include <vector>
#include <eigen-3.3.9/Eigen/Dense>
#include "simulationClasses.h"

using namespace std;
using namespace Eigen;

bool isContained(VectorXd& x0, vector<polytope>& polytopes){
    const double abs_tol = 1e-7;    // used by MPT to check containments.
    MatrixXd Aineq;
    VectorXd bineq;

    for (int i = 0; i < polytopes.size(); i++){
        Aineq = polytopes[i].A;
        bineq = polytopes[i].b;

    //     bool guard = true;
    //     int k = 0;
    //     while(guard && k < Aineq.rows()){
    //         if (Aineq.row(k)*x0 > bineq(k))
    //             if (abs(Aineq.row(k)*x0-bineq(k)) > abs_tol)
    //                 guard = false;
    //         k++;
    //     }

        VectorXd sat = Aineq * x0 - bineq;
        bool guard = (sat.array() <= abs_tol).all();
        if (guard)
            return true;
    }
    return false;
}


// Remove redundant inequalities of box constraints.
polytope simplify2box(VectorXd& x0, polytope& CIS, model& mdl, polytope& inputSet){
    const double abs_tol = 1e-7;   // used by MPT to check containments.
    // Get variables:
    MatrixXd Ad = mdl.Ad;
    MatrixXd Bd = mdl.Bd;
    MatrixXd cisA = CIS.A;
    VectorXd cisb = CIS.b;
    MatrixXd Gu = inputSet.A;
    MatrixXd Fu = inputSet.b;
    int ulen = Gu.cols();
    
    // Construct linear inequality constraints:
    // Aineq = cisA * Bd, bineq = cisb - cisA * Ad *x0
    MatrixXd Aineq(cisA.rows()+Gu.rows(), Bd.cols());
    Aineq << cisA * Bd, Gu;
    VectorXd bineq(cisb.rows()+Fu.rows());
    bineq << cisb - cisA * Ad * x0, Fu.col(0);

    // Extract the box:
    // box = [lb,ub]
    VectorXd ub = __DBL_MAX__ * VectorXd::Ones(Aineq.cols(),1);
    VectorXd lb = -__DBL_MAX__ * VectorXd::Ones(Aineq.cols(),1);
    // Construct box polytope:
    // box = { x | boxA x < boxb }, boxA = [-I; I], boxb = [-lb; ub]
    MatrixXd boxA(ulen*2, ulen);
    boxA << -MatrixXd::Identity(ulen,ulen), MatrixXd::Identity(ulen,ulen);
    MatrixXd boxb(ulen*2,1);

    // Normalize by the non-zero value of each row:
    bool guard = false;
    for (int j=0; j<Aineq.rows(); j++){
        int sum = 0;
        int idx;
        double val;
        for (int i=0; i<Aineq.cols(); i++){
            if (Aineq(j,i)!=0){
                val = Aineq(j,i);
                sum++;
                idx = i;
            }
        }
        if (sum>1){
            cout << "More than one non-zero value at row: "<< j << ", polytope not a hyper-rectangle" << endl;
            exit(-1);
        }
        else if (sum==1){
            val = abs(Aineq(j,idx));
            bineq(j) = bineq(j)/val;
            if (Aineq(j,idx)>0)
                ub(idx) = min(ub(idx),bineq(j));
            if (Aineq(j,idx)<0)
                lb(idx) = max(lb(idx),-bineq(j));
        }
        else    // Everything is zero, we check feasibility:
            if (bineq(j)<0 && abs(bineq(j)) > abs_tol) {
                // Infeasible: 
                ub = -VectorXd::Ones(Aineq.cols(),1);
                lb = VectorXd::Ones(Aineq.cols(),1);
                boxb << -lb, ub;
                polytope box(boxA,boxb);
                return box;
            }
    }
    
    // Add lb, ub to box and return:
    boxb << -lb, ub;
    polytope box(boxA,boxb);
    return box;
}

pair<double, VectorXd> callSupervisor(VectorXd& x_curr, VectorXd& u_des, model& mdl, polytope& RCIS, polytope& inputSet, int& method){
    // Store candidate solutions:
    double f_cand;
    VectorXd u_cand(u_des.rows());
    
    // RCIS inequalitites:
    MatrixXd rcisA = RCIS.A;
    VectorXd rcisb = RCIS.b;
    // wrt to u: rcisA * Bd * u < rcisb - rcisA * Ad * x_curr
    // Simplify to box constraints:
    polytope box = simplify2box(x_curr, RCIS, mdl, inputSet);
    MatrixXd Aineq = box.A;
    VectorXd bineq = box.b.col(0);

    if (method==1){   
        // Approach 1: Use Gurobi:
        // Original cost: || u - udes ||^2.
        MatrixXd H = MatrixXd::Identity(mdl.Nu,mdl.Nu);
        VectorXd c = -2 * u_des;

        // Call solver:
        // cout << "RCIS: " << k << endl;
        opt_result res = solveGurobi(H,c,Aineq,bineq,0);
        if (res.solved){
            u_cand = res.sol;
            f_cand = res.objVal;
        }
        else
            f_cand = __DBL_MAX__;
    }
    else{
        // Approach 2: Analytical solution:
        VectorXd lb(u_des.rows()), ub(u_des.rows());
        lb(0) = -bineq(0);   
        lb(1) = -bineq(1);   
        lb(2) = -bineq(2);
        ub(0) = bineq(3);   
        ub(1) = bineq(4);   
        ub(2) = bineq(5);

        if ((ub-lb).minCoeff()<0)   // infeasible
            f_cand = __DBL_MAX__;
        else{
            //           lb(j), if x0(j) < lb(j),
            // x*(k) =   ub(j), if x0(j) > ub(j),
            //           x0(j), otherwise.
            for (int j=0; j<ub.rows(); j++){
                if (u_des(j)<lb(j))
                    u_cand(j) = lb(j);
                else if (u_des(j)>ub(j))
                    u_cand(j) = ub(j);
                else
                    u_cand(j) = u_des(j);
            }
            f_cand = (u_cand-u_des).squaredNorm();
        }
    }
    
    pair<double, VectorXd> res(f_cand,u_cand);
    return res;
};
#endif
