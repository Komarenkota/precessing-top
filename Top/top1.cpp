#include <iostream>
#include <cmath>
#include <fstream>

double g = 9.81;

//std::ofstream outz("Fz.txt");


void vectorCrossProduct(double a[3],double b[3],double result[3]) {    
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

double computeTiltAngle(const double axis_world[3]) {
    double z_axis[3] = {0.0, 0.0, 1.0};
    
    double dot = axis_world[0]*z_axis[0] + axis_world[1]*z_axis[1] + axis_world[2]*z_axis[2];
    
    double norm_axis = sqrt(axis_world[0]*axis_world[0] + axis_world[1]*axis_world[1] + axis_world[2]*axis_world[2]);
    double norm_z = 1.0; 
    
    double angle = acos(dot / (norm_axis * norm_z));
    return angle;
}

void transpose(const double input[3][3], double output[3][3]) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            output[j][i] = input[i][j];  
        }
    }
}

void makeSkew(double input[3], double result[3][3]) {
    result[0][0] = 0;       result[0][1] = -input[2]; result[0][2] = input[1];
    result[1][0] = input[2]; result[1][1] = 0;        result[1][2] = -input[0];
    result[2][0] = -input[1]; result[2][1] = input[0]; result[2][2] = 0;
}

void addVectors(const double a[3], const double b[3], double result[3]) {
    for (int i = 0; i < 3; ++i) {
        result[i] = a[i] + b[i];    
    }
}

void subVectors(const double a[3], const double b[3], double result[3]) {
    for (int i = 0; i < 3; ++i) {
        result[i] = a[i]- b[i];    
    }
}

void subMatrices(const double a[3][3], const double b[3][3], double result[3][3]) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[i][j] = a[i][j]- b[i][j];
        }
    }
}

void addMatrices(const double a[3][3], const double b[3][3], double result[3][3]) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[i][j] = a[i][j]+ b[i][j];
        }
    }
}

void multiplyMatrices(const double a[3][3], const double b[3][3], double result[3][3]) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[i][j] = 0;
            for (int k = 0; k < 3; ++k) {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

void multiplyMatrixByVector(const double mat[3][3], const double vec[3], double result[3]) {
    for (int i = 0; i < 3; ++i) {
        result[i] = 0;
        for (int j = 0; j < 3; ++j) {
            result[i] += mat[i][j] * vec[j];
        }
    }
}

double scalar_multiplication(double* a, double* b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


void computeInertiaTensor(const double inv_I0[3][3], double R[3][3], double R_T[3][3], double inv_I[3][3] ){
    double temp2[3][3];
    multiplyMatrices(R,inv_I0, temp2);//R*inv(I0)
    multiplyMatrices(temp2, R_T, inv_I);//R*inv(I0)*R_T
}

typedef struct  // Структура кватерниона
{
    double s;  
    double v[3]; 
} quaternion;

void quaternion_multiplication(quaternion a, quaternion b, quaternion* res)
{
    res->s = a.s*b.s - scalar_multiplication(a.v, b.v);

    res->v[0] = a.s*b.v[0] + b.s*a.v[0] + a.v[1]*b.v[2] - a.v[2]*b.v[1];
    res->v[1] = a.s*b.v[1] + b.s*a.v[1] + a.v[2]*b.v[0] - a.v[0]*b.v[2];
    res->v[2] = a.s*b.v[2] + b.s*a.v[2] + a.v[0]*b.v[1] - a.v[1]*b.v[0];
}

void quaternion_normalize(quaternion* q)
{
    double l = sqrt(q->v[0]*q->v[0] + q->v[1]*q->v[1] + q->v[2]*q->v[2] + q->s*q->s);

    q->s = q->s / l;
    q->v[0] = q->v[0] / l;
    q->v[1] = q->v[1] / l;
    q->v[2] = q->v[2] / l; 
}

void quaternionToMatrix(quaternion q, double res[3][3])
{
    res[0][0] = 1 - 2*q.v[1]*q.v[1] - 2*q.v[2]*q.v[2];
    res[0][1] = 2*q.v[0]*q.v[1] - 2*q.s*q.v[2];
    res[0][2] = 2*q.v[0]*q.v[2] + 2*q.s*q.v[1];

    res[1][0] = 2*q.v[0]*q.v[1] + 2*q.s*q.v[2];
    res[1][1] = 1-2*q.v[0]*q.v[0]-2*q.v[2]*q.v[2];
    res[1][2] = 2*q.v[1]*q.v[2] - 2*q.s*q.v[0];
    
    res[2][0] = 2*q.v[0]*q.v[2] - 2*q.s*q.v[1];
    res[2][1] = 2*q.v[1]*q.v[2] + 2*q.s*q.v[0];
    res[2][2] = 1-2*q.v[0]*q.v[0]-2*q.v[1]*q.v[1];
}

struct Top
{
    const double m;            //mass
    const double I0[3][3];     //inertia tensor in body-space
    const double inv_I0[3][3]; // inverse of I0
    const double r0[3];        // r``

    double x[3];               //x(t)
    quaternion q;               /* q(t) */
    double P[3];               // P(t)
    double L[3];               // L(t)

    double inv_I[3][3];        // inverse of tensor inertia in world-space
    double R[3][3];                /* R(t) */
    double w[3];               // w(t)
    double v[3];               // v(t)

    double F[3];               // F(t)
    double torque[3];          // torque(t)

    double r[3];               // r` = R*r``
};


void checkEnergy(Top* top) {
    double E_rot = 0.5 * (top->w[0] * top->L[0] + top->w[1] * top->L[1] + top->w[2] * top->L[2]);
    printf("E =  %.15e\n", E_rot);
}

void StateToArray(Top *top, double *y)
{
    *y++ = top->x[0]; /* x component of position */
    *y++ = top->x[1]; 
    *y++ = top->x[2];
    
    *y++ = top->q.s;
    *y++ = top->q.v[0];
    *y++ = top->q.v[1];
    *y++ = top->q.v[2];

    *y++ = top->P[0];
    *y++ = top->P[1];
    *y++ = top->P[2];
    *y++ = top->L[0];
    *y++ = top->L[1];
    *y++ = top->L[2];
}


void ArrayToState(Top *top, double *y)
{
    top->x[0] = *y++;
    top->x[1] = *y++;
    top->x[2] = *y++;
    
    top->q.s = *y++;
    top->q.v[0] = *y++;
    top->q.v[1] = *y++;
    top->q.v[2] = *y++;

    top->P[0] = *y++;
    top->P[1] = *y++;
    top->P[2] = *y++;

    top->L[0] = *y++;
    top->L[1] = *y++;
    top->L[2] = *y++;

    top->v[0] = top->P[0] / top->m; // v(t) = P(t)/mass
    top->v[1] = top->P[1] / top->m;
    top->v[2] = top->P[2] / top->m;

    quaternion_normalize(&top->q);  // Нормализуем кватернион
    quaternionToMatrix(top->q, top->R);  // Считаем матрицу поворота

    double R_T[3][3];
    transpose(top->R, R_T);
    computeInertiaTensor(top->inv_I0, top->R, R_T, top->inv_I); // inv_I(t) = R(t)*inv_I0(t)*R_T(t) 
    multiplyMatrixByVector(top->inv_I, top->L, top->w); // w(t) = inv_I(t) * L(t)
}



double computeReactionForce(Top *top) {   
    double skew_w[3][3];
    makeSkew(top->w, skew_w);

    //double temp1[3][3];
    //multiplyMatrices(skew_w,top->inv_I,temp1);
    //double temp2[3][3];
    //transpose(skew_w, temp2);
    //double temp3[3][3];
    //multiplyMatrices(top->inv_I,temp2, temp3);
    //double temp4[3][3];
    //addMatrices(temp1,temp3,temp4);             /// w*inv(I) + inv(I)w*T
    //double temp5[3];
    //multiplyMatrixByVector(temp4,top->L, temp5);/// (w*inv(I) + inv(I)w*T)L
    double temp1 [3];
    vectorCrossProduct(top->L, top->w, temp1);
    double temp5[3];
    multiplyMatrixByVector(top->inv_I, temp1, temp5);   ///  inv(I)(L x w)

    
    double temp6[3][3];
    makeSkew(temp5, temp6);
    double temp7[3][3];
    multiplyMatrices(temp6,top->R, temp7);/// ((w*inv(I) + inv(I)w*T)L)*R
    double temp8[3];
    multiplyMatrixByVector(temp7, top->r0,temp8);   /// ((w*inv(I) + inv(I)w*T)L)*Rr``

    double temp9[3][3];
    multiplyMatrices(skew_w,top->R, temp9);
    double temp10[3][3];
    multiplyMatrices(skew_w,temp9, temp10);
    double temp11[3];
    multiplyMatrixByVector(temp10, top->r0, temp11); //w*w*Rr`` 
    
    double temp12[3];
    addVectors(temp11,temp8, temp12);               

    double numerator = g - temp12[2];
    double temp13[3];
    multiplyMatrixByVector(top->R,top->r0,temp13);      // r` = Rr``temp13[2] = 0;                                     
    double vec[3];
    double z[3] = {0.0,0.0,1.0};
    vectorCrossProduct(temp13, z, vec);                 // (r` x  z)
    double temp14[3];
    multiplyMatrixByVector(top->inv_I, vec, temp14);    // inv(I)(r` x  z)
    double temp15[3][3];
    makeSkew(temp14, temp15);                           // (inv(I)(r` x  z))*
    double temp16[3][3];
    multiplyMatrices(temp15, top->R, temp16);           // (inv(I)(r` x  z))*R
    double temp17[3];
    multiplyMatrixByVector(temp16, top->r0, temp17);    // (inv(I)(r` x  z))*Rr``


    double denominator = (1.0/top->m) + temp17[2];

    double Nz = numerator/denominator;
    return Nz;
}


void ComputeAccelerationOfContactPoint(Top *top, double Nz){
    double skew_w[3][3];
    makeSkew(top->w, skew_w);
    double N[3] = {0,0,Nz};

    // w*inv(I) + inv(I)w*T
    //`double temp1[3][3];
    //`multiplyMatrices(skew_w,top->inv_I,temp1);
    //`double temp2[3][3];
    //`transpose(skew_w, temp2);
    //`double temp3[3][3];
    //`multiplyMatrices(top->inv_I,temp2, temp3);
    //`double temp4[3][3];
    //`addMatrices(temp1,temp3,temp4);
    //`double temp5[3];
    //`multiplyMatrixByVector(temp4,top->L, temp5);/// (w*inv(I) + inv(I)w*T)L

    double temp1 [3];
    vectorCrossProduct(top->L, top->w, temp1);
    double temp5[3];
    multiplyMatrixByVector(top->inv_I, temp1, temp5);   ///  inv(I)(L x w)

    double temp6[3][3];
    makeSkew(temp5, temp6);
    double temp7[3][3];
    multiplyMatrices(temp6,top->R, temp7);/// ((w*inv(I) + inv(I)w*T)L)*R
    double temp8[3];
    multiplyMatrixByVector(temp7, top->r0,temp8);/// ((w*inv(I) + inv(I)w*)L)*Rr``


    double temp9[3][3];
    multiplyMatrices(skew_w,top->R, temp9);
    double temp10[3][3];
    multiplyMatrices(skew_w,temp9, temp10);
    double temp11[3];
    multiplyMatrixByVector(temp10, top->r0, temp11); //w*w*Rr`` 

    double temp12[3];
    multiplyMatrixByVector(top->R,top->r0,temp12);      // r` = Rr``
    double temp13[3];
    vectorCrossProduct(temp12, N, temp13);
    double temp14[3];
    multiplyMatrixByVector(top->inv_I,temp13, temp14); // inv(I) * (r` x N)
    double temp15[3][3];
    makeSkew(temp14, temp15);                       // (inv(I) * (r` x N))*
    double temp16[3][3];
    multiplyMatrices(temp15, top->R, temp16);
    double temp17[3];
    multiplyMatrixByVector(temp16, top->r0, temp17);

    double temp18[3];
    addVectors(temp11, temp8, temp18);
    addVectors(temp18, temp17, temp18);

    double temp19 = (N[2] / top->m) - g + temp18[2];

    std::cout <<"Ускорение точки опоры = " <<temp19<< '\n';
}

bool checkCollision(Top *top) {
    double world_r[3];
    double r[3];
    multiplyMatrixByVector(top->R, top->r0, r);
    addVectors(top->x, r, world_r);
    if (world_r[2] <= 0) {
 //       top->x[2] -= world_r[2]; 
        return true;
    }
    return false;
}

void computeTorqueAndForces(Top *top)
{
    double react_force[3] = {0.0,0.0,0.0};
    react_force[2] = computeReactionForce(top);
    std::cout << "Nz = "<<react_force[2]<<'\n';


    //outz << react_force[2] << '\n';

    ComputeAccelerationOfContactPoint(top, react_force[2]);

    double gravity_force[3];
    gravity_force[0] = 0;
    gravity_force[1] = 0;
    gravity_force[2] = -top->m*g;

    addVectors(react_force, gravity_force, top->F);

    //torque(t) = r`xFi  // r` = R(t)*r``
    double r[3];
    multiplyMatrixByVector(top->R,top->r0,r);
    vectorCrossProduct(r, react_force, top->torque);
}
        


void initState(Top* top) {
    double theta = 3.14 / 9.0;
    double x[3] = {0.0, 0.0, 0.1};
    quaternion q;
    q.s = cos(theta / 2.0);
    q.v[0] = sin(theta / 2.0);
    q.v[1] = 0.00;
    q.v[2] = 0.0;

    double P[3] = {0.0, 0.0, 0.0};
    double L[3] = {0.0, 0.0, 20.0};
    double y[13] = {
        x[0], x[1], x[2],      // Позиция
        q.s, q.v[0], q.v[1], q.v[2], // Кватернион
        P[0], P[1], P[2],      // Импульс
        L[0], L[1], L[2]       // Момент импульса
    };
    ArrayToState(top, y);

    // коррекция положения центра масс
    double offset[3];
    multiplyMatrixByVector(top->R, top->r0, offset);
    top->x[2] = -offset[2];

    StateToArray(top, y);
}

void DdtStateToArray(Top* tmp_top, double *current_state, double *xdot)
{   

/* copy d/dt x(t) = v(t) into xdot */
    *xdot++ = tmp_top->v[0];
    *xdot++ = tmp_top->v[1];
    *xdot++ = tmp_top->v[2];
   
    quaternion qdot;

    qdot.s = 0;
    qdot.v[0] = 0;
    qdot.v[1] = 0;
    qdot.v[2] = 0;


    quaternion omega_quat;

    omega_quat.s = 0;

    omega_quat.v[0] = tmp_top->w[0];
    omega_quat.v[1] = tmp_top->w[1];
    omega_quat.v[2] = tmp_top->w[2];

    quaternion_multiplication(omega_quat,tmp_top->q,  &qdot);   // dq/dt = 0.5*omega*q

    qdot.s *= 0.5;
    qdot.v[0] *= 0.5;
    qdot.v[1] *= 0.5;
    qdot.v[2] *= 0.5;
    
    *xdot++ = qdot.s;
    *xdot++ = qdot.v[0];
    *xdot++ = qdot.v[1];
    *xdot++ = qdot.v[2];

    *xdot++ = tmp_top->F[0]; /* d/dt P(t) = F(t) */
    *xdot++ = tmp_top->F[1];
    *xdot++ = tmp_top->F[2];
    
    *xdot++ = tmp_top->torque[0]; /* d/dt L(t) = τ(t) */
    *xdot++ = tmp_top->torque[1];
    *xdot++ = tmp_top->torque[2];
}




void rk4(Top* top, double* x0, double dt, double* x_fin) {
    const int n = 13;
    double k1[n], k2[n], k3[n], k4[n], temp[n];

    // K1
    {
        Top tmp = *top;
        ArrayToState(&tmp, x0);
        computeTorqueAndForces(&tmp);
        DdtStateToArray(&tmp, x0, k1);
    }

    // K2
    for (int i = 0; i < n; i++) {
        temp[i] = x0[i] + 0.5 * dt * k1[i];
    }
    {
        Top tmp = *top;
        ArrayToState(&tmp, temp);
        computeTorqueAndForces(&tmp);
        DdtStateToArray(&tmp, temp, k2);
    }

    // K3
    for (int i = 0; i < n; i++) {
        temp[i] = x0[i] + 0.5 * dt * k2[i];
    }
    {
        Top tmp = *top;
        ArrayToState(&tmp, temp);
        computeTorqueAndForces(&tmp);
        DdtStateToArray(&tmp, temp, k3);
    }

    // K4
    for (int i = 0; i < n; i++) {
        temp[i] = x0[i] + dt * k3[i];
    }
    {
        Top tmp = *top;
        ArrayToState(&tmp, temp);
        computeTorqueAndForces(&tmp);
        DdtStateToArray(&tmp, temp, k4);
    }

    for (int i = 0; i < n; i++) {
        x_fin[i] = x0[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }

    ArrayToState(top, x_fin);
}