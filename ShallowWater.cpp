#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "ShallowWater.h"

#define SW ShallowWater
#define f1(u,v,du_dx,du_dy,dh_dx) (-(u*du_dx + v*du_dy + g*dh_dx ))
#define f2(u,v,dv_dx,dv_dy,dh_dy) (-(u*dv_dx + v*dv_dy + g*dh_dy ))
#define f3(u,v,h,du_dx,dh_dx,dv_dy,dh_dy) (-(u*dh_dx + h*du_dx + v*dh_dy + h*dv_dy))

#define F77NAME(x) x##_

extern "C"{
    
    void F77NAME(dgemv) (const char &trans, const int &m, const int &n, const double &alpha, const double *a, 
    const int &lda, const double *x, const int &incx, const double &beta, double *y, const int &incy);
    
    void F77NAME(dgemm) (const char &trans, const char &transb, const int &m1, const int &n, const int &k, const double &alpha,
    const double *a, const int &lda, const double *b, const int &ldb, const double &beta, double *c, const int &ldc);
    
    void F77NAME(dcopy) (const int &n, const double *x, const int &incx, double *y, const int &incy);
    
    void F77NAME(daxpy) (const int &n, const double &alpha, const double *x, const int &incx, const double *y, const int &incy);
}

SW::ShallowWater  (double dt, int T, int Nx, int Ny, 
                        int ic, double dx, double dy) {

    this->dt   = dt;
    this->T    = T;
    this->Nx   = Nx;
    this->Ny   = Ny;
    this->ic   = ic;
    this->dx   = dx;
    this->dy   = dy;

    // Allocate memory for the matrix, initialise to zero
    try {
        U1 = new double[Nx*Ny]();
        U2 = new double[Nx*Ny]();
        V1 = new double[Nx*Ny]();
        V2 = new double[Nx*Ny]();
        H1 = new double[Nx*Ny]();
        H2 = new double[Nx*Ny]();
        du_dx = new double[Nx*Ny]();
        du_dy = new double[Nx*Ny]();
        dv_dx = new double[Nx*Ny]();
        dv_dy = new double[Nx*Ny]();
        dh_dx = new double[Nx*Ny]();
        dh_dy = new double[Nx*Ny]();
        dU_dt = new double[Nx*Ny]();
        dV_dt = new double[Nx*Ny]();
        dH_dt = new double[Nx*Ny]();
        k2_U = new double[Nx*Ny]();
        k2_V = new double[Nx*Ny]();
        k2_H = new double[Nx*Ny]();
        k3_U = new double[Nx*Ny]();
        k3_V = new double[Nx*Ny]();
        k3_H = new double[Nx*Ny]();
        k4_U = new double[Nx*Ny]();
        k4_V = new double[Nx*Ny]();
        k4_H = new double[Nx*Ny]();
        CoeffX = new double[Nx*Nx]();
        CoeffY = new double[Ny*Ny]();
    } catch (std::bad_alloc& ex) {
        std::cout << "Out of memory!" << std::endl;
        std::exit(1);
    }
}

SW::~ShallowWater() {
    if(U1) delete U1;
    if(U2) delete U2;
    if(V1) delete V1;
    if(V2) delete V2;
    if(H1) delete H1;
    if(H2) delete H2;
    if(du_dx) delete du_dx;
    if(du_dy) delete du_dy;
    if(dv_dx) delete dv_dx;
    if(dv_dy) delete dv_dy;
    if(dh_dx) delete dh_dx;
    if(dh_dy) delete dh_dy;
    if(dU_dt) delete dU_dt;
    if(dV_dt) delete dV_dt;
    if(dH_dt) delete dH_dt;
    if(k2_U ) delete k2_U;
    if(k2_V ) delete k2_V;
    if(k2_H ) delete k2_H;
    if(k3_U ) delete k3_U;
    if(k3_V ) delete k3_V;
    if(k3_H ) delete k3_H;
    if(k4_U ) delete k4_U;
    if(k4_V ) delete k4_V;
    if(k4_H ) delete k4_H;
    if(CoeffX) delete CoeffX;
    if(CoeffY) delete CoeffY;
    
}

void SW::SetInitialConditions() { 

    // Set Initial conditions for U and V
    /*for (int i = 0; i<Nx*Ny; ++i) {
        V1[i] = 0.0;
        U1[i] = 0.0;  //probably already initialised to 0 above when creating dynamic 
    }
    */
    switch (ic) {
        case 1: // test case 1
            for (int i = 0; i<(Nx); ++i) {
                double exp1 = ((dx*i)-50.0);
                for (int j = 0; j<Ny; ++j) {
                    
                    H1[j+Ny*i] = 10.0 + exp(-(exp1*exp1)*0.04);
                }           
            }
            break;
        case 2: // test case 2
            for (int i = 0; i<(Nx); ++i) {
                for (int j = 0; j<Ny; ++j) {
                    double exp2 = ((dy*j)-50.0);
                    H1[j+Ny*i] = 10.0 + exp(-(exp2*exp2)*0.04);
                }           
            }
            break;
        case 3: // test case 3
            for (int i = 0; i<(Nx); ++i) {
                for (int j = 0; j<Ny; ++j) {
                    H1[j+Ny*i] = 10.0 + exp(-((dx*i-50.0)*(dx*i-50.0) + (dy*j-50.0)*(dy*j-50.0))*0.04);
                }           
            }
            break;
        case 4: // test case 4
            for (int i = 0; i<(Nx); ++i) {
                for (int j = 0; j<Ny; ++j) {
                    H1[j+Ny*i] = 10.0 + exp(-((dx*i-25.0)*(dx*i-25.0) + (dy*j-25.0)*(dy*j-25.0))*0.04) + 
                                 exp(-((dx*i-75.0)*(dx*i-75.0) + (dy*j-75.0)*(dy*j-75.0))*0.04);
                }           
            }
            break;   
    }
}

void SW::TimeIntegrate(int option) {
    // Cannot parallelise this, due to race condition
    
    switch(option) {
        case 1:
            for (double t = dt; t<=T; t += dt) {
                TimeIntegrateSingle();
            }      
            break;
        
        case 2:
            SW::FillMatrix();
            for (double t = dt; t<=T; t += dt) {
                BLASTimeIntegrateSingle();
            }    
            break;
    }   
  
}

void SW::TimeIntegrateSingle() {

    // Copy to local scope, avoid referencing
    double ddt = dt;
    int Nxx = Nx;
    int Nyy = Ny;
    double const1 = -1.0/60.0;
    const double g = 9.81;
     
    SW::gradient(U1,V1, H1);
    
    double dt2 = ddt*0.5;
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < Nxx*Nyy; ++i) {
        // Evaluate k1
        dU_dt[i] = f1(U1[i], V1[i], du_dx[i], du_dy[i], dh_dx[i]); // Evaluate the time derivative of U using the shallow water equations
        dV_dt[i] = f2(U1[i], V1[i], dv_dx[i], dv_dy[i], dh_dy[i]); // Evaluate the time derivative of V using the shallow water equations
        dH_dt[i] = f3(U1[i], V1[i], H1[i], du_dx[i], dh_dx[i], dv_dy[i], dh_dy[i]); 
        
        U2[i] = U1[i] + dU_dt[i]*dt2;
        V2[i] = V1[i] + dV_dt[i]*dt2;
        H2[i] = H1[i] + dH_dt[i]*dt2;
    }
    
    SW::gradient(U2,V2, H2); // caclulate gradient with new U,V,H  // calculates new du/dx
    
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < Nxx*Nyy; ++i) {
    
        k2_U[i] = f1(U2[i], V2[i], du_dx[i], du_dy[i], dh_dx[i]); // Evaluate the time derivative of U using the shallow water equations
        k2_V[i] = f2(U2[i], V2[i], dv_dx[i], dv_dy[i], dh_dy[i]); // Evaluate the time derivative of V using the shallow water equations
        k2_H[i] = f3(U2[i], V2[i], H2[i], du_dx[i], dh_dx[i], dv_dy[i], dh_dy[i]); 

        U2[i] = U1[i] + k2_U[i]*dt2;
        V2[i] = V1[i] + k2_V[i]*dt2;
        H2[i] = H1[i] + k2_H[i]*dt2;
    }
    
    SW::gradient(U2, V2, H2); // caclulate gradient with new U,V,H  // calculates new du/dx
    
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < Nxx*Nyy; ++i) {
        
        k3_U[i] = f1(U2[i], V2[i], du_dx[i], du_dy[i], dh_dx[i]); // Evaluate the time derivative of U using the shallow water equations
        k3_V[i] = f2(U2[i], V2[i], dv_dx[i], dv_dy[i], dh_dy[i]); // Evaluate the time derivative of V using the shallow water equations
        k3_H[i] = f3(U2[i], V2[i], H2[i], du_dx[i], dh_dx[i], dv_dy[i], dh_dy[i]); 

        U2[i] = U1[i] + k3_U[i]*ddt;
        V2[i] = V1[i] + k3_V[i]*ddt;
        H2[i] = H1[i] + k3_H[i]*ddt;
        
    }
    
    SW::gradient(U2, V2,H2); // caclulate gradient with new U,V,H  // calculates new du/dx
    
    #pragma omp parallel for schedule(static) 
    for (int i = 0; i < Nxx*Nyy; ++i) {
        
        k4_U[i] = f1(U2[i], V2[i], du_dx[i], du_dy[i], dh_dx[i]); // Evaluate the time derivative of U using the shallow water equations
        k4_V[i] = f2(U2[i], V2[i], dv_dx[i], dv_dy[i], dh_dy[i]); // Evaluate the time derivative of V using the shallow water equations
        k4_H[i] = f3(U2[i], V2[i], H2[i], du_dx[i], dh_dx[i], dv_dy[i], dh_dy[i]); 
      
    }
    double constF = -10*const1;
    #pragma omp parallel for schedule(static) 
    for (int i = 0; i < Nxx*Nyy; ++i) {
        // Update the solution using RK4
        
        U1[i] += (constF)*(dU_dt[i] + 2*k2_U[i] + 2*k3_U[i] + k4_U[i])*ddt;
        V1[i] += (constF)*(dV_dt[i] + 2*k2_V[i] + 2*k3_V[i] + k4_V[i])*ddt;
        H1[i] += (constF)*(dH_dt[i] + 2*k2_H[i] + 2*k3_H[i] + k4_H[i])*ddt;
        
    }
}


void SW::gradient(double* U1, double* V1, double* H1){
    
    double ddx = dx;
    double ddy = dy;
    int Nxx = Nx;
    int Nyy = Ny;
    double invDx = 1.0/ddx;
    double invDy = 1.0/ddy;
    double const1 = -1.0/60.0;
    
    #pragma omp parallel for schedule(static) collapse(2)
    for (int i = 0; i < Nxx; ++i) {
        for (int j = 0; j < Nyy; ++j) { // this is column major order, so goes dowm column first
            int index = j+Nyy*i;
            
            
            du_dx[index] = invDx*(  const1* (i <     3 ? U1[index+(Nxx-3)*Nyy] : U1[index   - 3*Nyy]) + //
                                    0.15  * (i <     2 ? U1[index+(Nxx-2)*Nyy] : U1[index   - 2*Nyy]) - // if U_i-n, the i-n < 0 where n = 1,2,3
                                    0.75  * (i <     1 ? U1[index+(Nxx-1)*Nyy] : U1[index   -   Nyy]) + //
                                    0.75  * (i > (Nxx-1)-1 ? U1[index-(Nxx-  1)*Nyy] : U1[index +   Nyy]) - //
                                    0.15  * (i > Nxx-2-1 ? U1[index-(Nxx-  2)*Nyy] : U1[index + 2*Nyy]) - // if U_i+n, the i+n > Nxx where n = 1,2,3
                                    const1* (i > Nxx-3-1 ? U1[index-(Nxx-  3)*Nyy] : U1[index + 3*Nyy])   //
                                    );
            // do dv/dx
            dv_dx[index] = invDx*(  const1* (i <     3 ? V1[index+(Nxx-3)*Nyy] : V1[index   - 3*Nyy]) + //
                                    0.15  * (i <     2 ? V1[index+(Nxx-2)*Nyy] : V1[index   - 2*Nyy]) - // if U_i-n, the i-n < 0 where n = 1,2,3
                                    0.75  * (i <     1 ? V1[index+(Nxx-1)*Nyy] : V1[index   -   Nyy]) + //
                                    0.75  * (i > (Nxx-1)-1 ? V1[index-(Nxx-  1)*Nyy] : V1[index +   Nyy]) - //
                                    0.15  * (i > Nxx-2-1 ? V1[index-(Nxx-  2)*Nyy] : V1[index + 2*Nyy]) - // if U_i+n, the i+n > Nxx where n = 1,2,3
                                    const1* (i > Nxx-3-1 ? V1[index-(Nxx-  3)*Nyy] : V1[index + 3*Nyy])   //
                                    );
            // do dh/dx
            dh_dx[index] = invDx*(  const1* (i <     3 ? H1[index+(Nxx-3)*Nyy] : H1[index   - 3*Nyy]) + //
                                    0.15  * (i <     2 ? H1[index+(Nxx-2)*Nyy] : H1[index   - 2*Nyy]) - // if U_i-n, the i-n < 0 where n = 1,2,3
                                    0.75  * (i <     1 ? H1[index+(Nxx-1)*Nyy] : H1[index   -   Nyy]) + //
                                    0.75  * (i > (Nxx-1)-1 ? H1[index-(Nxx-  1)*Nyy] : H1[index +   Nyy]) - //
                                    0.15  * (i > Nxx-2-1 ? H1[index-(Nxx-  2)*Nyy] : H1[index + 2*Nyy]) - // if U_i+n, the i+n > Nxx where n = 1,2,3
                                    const1* (i > Nxx-3-1 ? H1[index-(Nxx-  3)*Nyy] : H1[index + 3*Nyy])   //
                                    );
            
            du_dy[index] = invDy*(  const1* (j <     3 ? U1[index+(Nyy-3)] : U1[index - 3]) + //
                                    0.15  * (j <     2 ? U1[index+(Nyy-2)] : U1[index - 2]) - // if U_i-n, the i-n < 0 where n = 1,2,3
                                    0.75  * (j <     1 ? U1[index+(Nyy-1)] : U1[index - 1]) + //
                                    0.75  * (j > (Nyy-1)-1 ? U1[index-Nyy+  1] : U1[index + 1]) - //
                                    0.15  * (j > Nyy-2-1 ? U1[index-Nyy+  2] : U1[index + 2]) - // if U_i+n, the i+n > Nxx where n = 1,2,3
                                    const1* (j > Nyy-3-1 ? U1[index-Nyy+  3] : U1[index + 3])   //
                                    );
                                    
            dv_dy[index] = invDy*(  const1* (j <     3 ? V1[index+(Nyy-3)] : V1[index - 3]) + //
                                    0.15  * (j <     2 ? V1[index+(Nyy-2)] : V1[index - 2]) - // if U_i-n, the i-n < 0 where n = 1,2,3
                                    0.75  * (j <     1 ? V1[index+(Nyy-1)] : V1[index - 1]) + //
                                    0.75  * (j > (Nyy-1)-1 ? V1[index-Nyy+  1] : V1[index + 1]) - //
                                    0.15  * (j > Nyy-2-1 ? V1[index-Nyy+  2] : V1[index + 2]) - // if U_i+n, the i+n > Nxx where n = 1,2,3
                                    const1* (j > Nyy-3-1 ? V1[index-Nyy+  3] : V1[index + 3])   //
                                    );
            
            dh_dy[index] = invDy*(  const1* (j <     3 ? H1[index+(Nyy-3)] : H1[index - 3]) + //
                                    0.15  * (j <     2 ? H1[index+(Nyy-2)] : H1[index - 2]) - // if U_i-n, the i-n < 0 where n = 1,2,3
                                    0.75  * (j <     1 ? H1[index+(Nyy-1)] : H1[index - 1]) + //
                                    0.75  * (j > (Nyy-1)-1 ? H1[index-Nyy+  1] : H1[index + 1]) - //
                                    0.15  * (j > Nyy-2-1 ? H1[index-Nyy+  2] : H1[index + 2]) - // if U_i+n, the i+n > Nxx where n = 1,2,3
                                    const1* (j > Nyy-3-1 ? H1[index-Nyy+  3] : H1[index + 3])   //
                                    );
        }
    }
    
}

void SW::FillMatrix() {
    double ddx = dx;
    // double ddy = dy;
    int Nyy = Ny;
    int Nxx = Nx;
    double const1 = 1.0/60.0/ddx;
    double const2 =     0.15/ddx;
    double const3 =     0.75/ddx;
    
    
    /*for (int i=0; i<Nxx*Nxx; i++) {
        CoeffX[i] = 0;
    }*/
    
    CoeffX[Nxx*(Nxx-1)  ] = -const3; // top right
    CoeffX[Nxx*(Nxx-2)  ] =  const2; // top right then left 
    CoeffX[Nxx*(Nxx-3)  ] = -const1; // top right then left then left
    CoeffX[Nxx*(Nxx-1)+1] =  const2; // top right then down
    CoeffX[Nxx*(Nxx-1)+2] = -const1; // top right then down then down
    CoeffX[Nxx*(Nxx-2)+1] = -const1; // second top right
    
    CoeffX[Nxx-1] = const3; // bottom left
    CoeffX[Nxx-2] =  -const2; // bottom left up  
    CoeffX[Nxx-3] = const1; // bottom left up then up
    CoeffX[2*Nxx-1] = -const2; // bottom left then right 
    CoeffX[2*Nxx-2] = const1; // bottom left then right then up
    CoeffX[3*Nxx-1] = const1; // bottom left then right then right
    
    // fills the 0.75 and -0.75
     
    for (int i = 1; i<Nxx; ++i) {
        int index = i*Nxx+i; // this index goes along the diagonal
        CoeffX[index-1  ] =  const3;
        CoeffX[index-Nxx] = -const3;
    }
    
    
    for (int i = 2; i<Nxx; ++i) {
        int index = i*Nxx+i; // this index goes along the diagonal
        CoeffX[index-2  ] = -const2;
        CoeffX[index-2*Nxx] =  const2;
    }
    
    
    for (int i = 3; i<Nxx; ++i) {
        int index = i*Nxx+i; // this index goes along the diagonal
        CoeffX[index-3  ] = const1;
        CoeffX[index-3*Nxx] =  -const1;
    }
    
    // now need matrix for vertical scenario
    
    for (int i=0; i<Nyy*Nyy; ++i) {
        CoeffY[i] = 0;
    }
            
    // this pre defined should be opposite sign to the x case 
    CoeffY[Nyy*(Nyy-1)  ] = -const3; // top right
    CoeffY[Nyy*(Nyy-2)  ] = const2; // top right then left 
    CoeffY[Nyy*(Nyy-3)  ] = -const1; // top right then left then left
    CoeffY[Nyy*(Nyy-1)+1] =  const2; // top right then down
    CoeffY[Nyy*(Nyy-1)+2] = -const1; // top right then down then down
    CoeffY[Nyy*(Nyy-2)+1] = -const1; // second top right
    
    CoeffY[Nyy-1] = const3; // bottom left
    CoeffY[Nyy-2] =  -const2; // bottom left up  
    CoeffY[Nyy-3] = const1; // bottom left up then up
    CoeffY[2*Nyy-1] = -const2; // bottom left then right 
    CoeffY[2*Nyy-2] = const1; // bottom left then right then up
    CoeffY[3*Nyy-1] = const1; // bottom left then right then right
    
    // fills the 0.75 and -0.75
    
    for (int i = 1; i<Nyy; ++i) {
        int index = i*Nyy+i; // this index goes along the diagonal
        CoeffY[index-1  ] =  const3;
        CoeffY[index-Nyy] = -const3;
    }
    
    
    for (int i = 2; i<Nyy; ++i) {
        int index = i*Nyy+i; // this index goes along the diagonal
        CoeffY[index-2  ] = -const2;
        CoeffY[index-2*Nyy] =  const2;
    }
    
    for (int i = 3; i<Nyy; ++i) {
        int index = i*Nyy+i; // this index goes along the diagonal
        CoeffY[index-3  ] = const1;
        CoeffY[index-3*Nyy] =  -const1;
    }
    
}

void SW::BLASGrad(double* U, double* V, double* H) {
    
    int Nxx = Nx;
    int Nyy = Ny;
    
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < Nyy; ++i){
        
        F77NAME(dgemv)('N', Nxx, Nxx, 1.0, CoeffX, Nxx, U + i, Nyy, 0.0, du_dx + i, Nyy);
        
        F77NAME(dgemv)('N', Nxx, Nxx, 1.0, CoeffX, Nxx, V + i, Nyy, 0.0, dv_dx + i, Nyy);
        
        F77NAME(dgemv)('N', Nxx, Nxx, 1.0, CoeffX, Nxx, H + i, Nyy, 0.0, dh_dx + i, Nyy);
    }
    
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < Nxx; ++j){
        F77NAME(dgemv)('N', Nyy, Nyy, 1.0, CoeffY, Nyy, U + (j*Nyy), 1, 0.0, du_dy + (j*Nyy), 1);
        
        F77NAME(dgemv)('N', Nyy, Nyy, 1.0, CoeffY, Nyy, V + (j*Nyy), 1, 0.0, dv_dy + (j*Nyy), 1);
        
        F77NAME(dgemv)('N', Nyy, Nyy, 1.0, CoeffY, Nyy, H + (j*Nyy), 1, 0.0, dh_dy + (j*Nyy), 1);
    }
}


void SW::BLASSolve(){
    int Nxx = Nx;
    int Nyy = Ny;
    double ddt = dt;
    const double g = 9.81;
    
    double U_val;
    double V_val;
    double H_val;
    
    // V2 is what will be used to add at the end yn + 1/6 (...) 
    
    double dudx_val;
    double dudy_val;
    double dvdx_val;
    double dvdy_val;
    double dhdx_val;
    double dhdy_val;
    
    
    
    SW::BLASGrad(U1,V1,H1);
    
    for (int i = 0; i < Nxx*Nyy; ++i) {
        // Evaluate k1
        U_val = U1[i];
        V_val = V1[i];
        H_val = H1[i];
        
        // V2 is what will be used to add at the end yn + 1/6 (...) 
        
        dudx_val = du_dx[i];
        dudy_val = du_dy[i];
        dvdx_val = dv_dx[i];
        dvdy_val = dv_dy[i];
        dhdx_val = dh_dx[i];
        dhdy_val = dh_dy[i];
        
        dU_dt[i] = f1(U_val,V_val,dudx_val,dudy_val,dhdx_val); // Evaluate the time derivative of U using the shallow water equations
        dV_dt[i] = f2(U_val,V_val,dvdx_val,dvdy_val,dhdy_val); // Evaluate the time derivative of V using the shallow water equations
        dH_dt[i] = f3(U_val,V_val,H_val,dudx_val,dhdx_val,dvdy_val,dhdy_val); // Evaluate the time derivative of H using the shallow w
        
    }
   
    F77NAME(dcopy) (Nxx*Nyy, U1, 1, U2, 1);     // copy U1 into U2
    F77NAME(daxpy) (Nxx*Nyy, ddt*0.5, dU_dt, 1, U2, 1); // U2 = U2 + k1*dt/2
  
    F77NAME(dcopy) (Nxx*Nyy, V1, 1, V2, 1);
    F77NAME(daxpy) (Nxx*Nyy, ddt*0.5, dV_dt, 1, V2, 1);

    F77NAME(dcopy) (Nxx*Nyy, H1, 1, H2, 1);
    F77NAME(daxpy) (Nxx*Nyy, ddt*0.5, dH_dt, 1, H2, 1);
    
    SW::BLASGrad(U2,V2,H2);
    
    for (int i = 0; i < Nxx*Nyy; ++i) {
        // Evaluate k2
        U_val = U2[i];
        V_val = V2[i];
        H_val = H2[i];
        
        // V2 is what will be used to add at the end yn + 1/6 (...) 
        
        dudx_val = du_dx[i];
        dudy_val = du_dy[i];
        dvdx_val = dv_dx[i];
        dvdy_val = dv_dy[i];
        dhdx_val = dh_dx[i];
        dhdy_val = dh_dy[i];
        
        k2_U[i] = f1(U_val,V_val,dudx_val,dudy_val,dhdx_val); // Evaluate the time derivative of U using the shallow water equations
        k2_V[i] = f2(U_val,V_val,dvdx_val,dvdy_val,dhdy_val); // Evaluate the time derivative of V using the shallow water equations
        k2_H[i] = f3(U_val,V_val,H_val,dudx_val,dhdx_val,dvdy_val,dhdy_val); // Evaluate the time derivative of H using the shallow w

    }
    
    F77NAME(dcopy) (Nxx*Nyy, U1, 1, U2, 1);
    F77NAME(daxpy) (Nxx*Nyy, ddt*0.5, k2_U, 1, U2, 1);
    F77NAME(dcopy) (Nxx*Nyy, V1, 1, V2, 1);
    F77NAME(daxpy) (Nxx*Nyy, ddt*0.5, k2_V, 1, V2, 1);       
    F77NAME(dcopy) (Nxx*Nyy, H1, 1, H2, 1);
    F77NAME(daxpy) (Nxx*Nyy, ddt*0.5, k2_H, 1, H2, 1);
    
    SW::BLASGrad(U2,V2,H2);

    for (int i = 0; i < Nxx*Nyy; ++i) {
        // Evaluate k3
        U_val = U2[i];
        V_val = V2[i];
        H_val = H2[i];
        
        dudx_val = du_dx[i];
        dudy_val = du_dy[i];
        dvdx_val = dv_dx[i];
        dvdy_val = dv_dy[i];
        dhdx_val = dh_dx[i];
        dhdy_val = dh_dy[i];
        
        k3_U[i] = f1(U_val,V_val,dudx_val,dudy_val,dhdx_val); // Evaluate the time derivative of U using the shallow water equations
        k3_V[i] = f2(U_val,V_val,dvdx_val,dvdy_val,dhdy_val); // Evaluate the time derivative of V using the shallow water equations
        k3_H[i] = f3(U_val,V_val,H_val,dudx_val,dhdx_val,dvdy_val,dhdy_val); // Evaluate the time derivative of H using the shallow w
        

    }
    
    F77NAME(dcopy) (Nxx*Nyy, U1,   1, U2, 1);
    F77NAME(daxpy) (Nxx*Nyy, ddt, k3_U, 1, U2, 1);
    
    F77NAME(dcopy) (Nxx*Nyy, V1,   1, V2, 1);
    F77NAME(daxpy) (Nxx*Nyy, ddt, k3_V, 1, V2, 1);
    
    F77NAME(dcopy) (Nxx*Nyy, H1,   1, H2, 1);
    F77NAME(daxpy) (Nxx*Nyy, ddt, k3_H, 1, H2, 1);
    
    SW::BLASGrad(U2,V2,H2);
    
    for (int i = 0; i < Nxx*Nyy; ++i) {
        // Evaluate k3
        U_val = U2[i];
        V_val = V2[i];
        H_val = H2[i];
        
        dudx_val = du_dx[i];
        dudy_val = du_dy[i];
        dvdx_val = dv_dx[i];
        dvdy_val = dv_dy[i];
        dhdx_val = dh_dx[i];
        dhdy_val = dh_dy[i];
        
        k4_U[i] = f1(U_val,V_val,dudx_val,dudy_val,dhdx_val); // Evaluate the time derivative of U using the shallow water equations
        k4_V[i] = f2(U_val,V_val,dvdx_val,dvdy_val,dhdy_val); // Evaluate the time derivative of V using the shallow water equations
        k4_H[i] = f3(U_val,V_val,H_val,dudx_val,dhdx_val,dvdy_val,dhdy_val); // Evaluate the time derivative of H using the shallow w
    }
    
}

void SW::BLASRK4() {
    //Precompute fractions
    double ddt = dt;
    double Nxx = Nx;
    double Nyy = Ny;
    double const1 = ddt/6.0;
    double const2 = const1*2.0;
   
    // u 
    F77NAME(daxpy) (Nxx*Nyy, const1, dU_dt, 1, U1, 1);  // k1 term in RK4
    F77NAME(daxpy) (Nxx*Nyy, const2, k2_U, 1, U1, 1); // k2 term in RK4
    F77NAME(daxpy) (Nxx*Nyy, const2, k3_U, 1, U1, 1); // k3 term in RK4
    F77NAME(daxpy) (Nxx*Nyy, const1, k4_U, 1, U1, 1);   // k4 term in RK4
    // v
    F77NAME(daxpy) (Nxx*Nyy, const1, dV_dt, 1, V1, 1);
    F77NAME(daxpy) (Nxx*Nyy, const2, k2_V, 1, V1, 1); 
    F77NAME(daxpy) (Nxx*Nyy, const2, k3_V, 1, V1, 1);
    F77NAME(daxpy) (Nxx*Nyy, const1, k4_V, 1, V1, 1); 
    // h
    F77NAME(daxpy) (Nxx*Nyy, const1, dH_dt, 1, H1, 1);
    F77NAME(daxpy) (Nxx*Nyy, const2, k2_H, 1, H1, 1); 
    F77NAME(daxpy) (Nxx*Nyy, const2, k3_H, 1, H1, 1);
    F77NAME(daxpy) (Nxx*Nyy, const1, k4_H, 1, H1, 1);     
}
    
void SW::BLASTimeIntegrateSingle() {
    
    SW::BLASSolve();
    
    SW::BLASRK4();
    
}

void SW::writeOutput() {
    std::ofstream outputFile("output.txt");
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            outputFile << std::setw(5) << j << std::setw(5) << i << std::setw(15) << U1[j*Ny+i] << std::setw(15) << V1[j*Ny+i] << std::setw(15) << H1[j*Ny+i] <<std::endl;
        }
        outputFile << std::endl;
    }
    outputFile.close();
}