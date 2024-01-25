#ifndef SHALLOWWATER_H
#define SHALLOWWATER_H

class ShallowWater {
    private:
        int Nx, Ny, T, ic;
        double dt,dx,dy;
        double* U1;
        double* V1;
        double* U2;
        double* V2;
        double* H1;
        double* H2;
        double* du_dx;
        double* du_dy;
        double* dv_dx;
        double* dv_dy;
        double* dh_dx;
        double* dh_dy;
        double* dU_dt;
        double* dV_dt;
        double* dH_dt;
        double* k2_U;
        double* k2_V;
        double* k2_H;
        double* k3_U;
        double* k3_V;
        double* k3_H;
        double* k4_U;
        double* k4_V;
        double* k4_H;
        double* CoeffX;
        double* CoeffY;
        /**
         * @brief Calculates the spatial gradients for the for loop method.
         * 
         * @param U     U matrix
         * @param V     V matrix
         * @param H1    H matrix 
         */
        void gradient(double* U1, double* V1, double* H1);
        
        /**
         * @brief computes the time derivatives thats needed for the k values in Runge-Kutta
         */
        void BLASSolve();
        /**
         * @brief Calculates the spatial gradients using BLAS
         * @param U1    U matrix
         * @param V1    V matrix
         * @param H1    H matrix
         */
        void BLASGrad(double* U1, double* V1, double* H1);
        /**
         * @brief Fills the stencil matrix where there are two scenarios for x and y
         */
        void FillMatrix();
        /**
         * @brief add the terms for the Runge-Kutta to calculate y_(n+1)
         */
        void BLASRK4();
        
    public:
        /**
         * @brief Set the Parameters required for the problem, initialise the 2D matrix.
         * 
         * @param dt    Time-step to use.
         * @param T     Total integration time.
         * @param Nx    Number of grid points in x
         * @param Ny    Number of grid points in y
         * @param ic    Index of the initial condition to use (1-4
         */
        ShallowWater(double dt, int T, int Nx, int Ny, int ic,
                          double dx = 1, double dy = 1);
        /// Destructor, free memory
        ~ShallowWater();

        /**
         * @brief Set the Initial Conditions of the problem
         */
        void SetInitialConditions();

        /**
         * @brief Perform the full solution.
         * @param option to choose For loop or BLAS method
         */
        void TimeIntegrate(int option);

        /**
         * @brief Perform one time step of the integration using for loop
         */
        void TimeIntegrateSingle();
        
        /**
         * @brief Perform one time step of the integration using BLAS
         */
        void BLASTimeIntegrateSingle();

        /**
         * @brief Output entire grid.
         */
        void writeOutput();
        
};

#endif