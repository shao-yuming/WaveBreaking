#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define DIMENSION 3

typedef struct
{
    double x, y, z;
} vec3;

typedef double newmat3[DIMENSION][DIMENSION];

void matrixMultiply(double A[DIMENSION][DIMENSION], double B[DIMENSION][DIMENSION], double result[DIMENSION][DIMENSION])
{
    for (int i = 0; i < DIMENSION; i++)
    {
        for (int j = 0; j < DIMENSION; j++)
        {
            result[i][j] = 0.0;
            for (int k = 0; k < DIMENSION; k++)
            {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void transpose(double matrix[DIMENSION][DIMENSION], double result[DIMENSION][DIMENSION])
{
    for (int i = 0; i < DIMENSION; i++)
    {
        for (int j = 0; j < DIMENSION; j++)
        {
            result[i][j] = matrix[j][i];
        }
    }
}



void getOmegaRTransformMatrix(vec3 unitEigenVector, newmat3 Q)
{
    vec3 r = unitEigenVector;
    vec3 vz = {0, 0, 1};

    vec3 gamma = {0, 0, 0};
    double cross_vz_r_x = vz.y * r.z - vz.z * r.y;
    double cross_vz_r_y = vz.z * r.x - vz.x * r.z;
    double cross_vz_r_z = vz.x * r.y - vz.y * r.x;
    double norm_cross_vz_r = sqrt(cross_vz_r_x * cross_vz_r_x + cross_vz_r_y * cross_vz_r_y + cross_vz_r_z * cross_vz_r_z);

    if (norm_cross_vz_r != 0)
    {
        gamma.x = cross_vz_r_x / norm_cross_vz_r;
        gamma.y = cross_vz_r_y / norm_cross_vz_r;
        gamma.z = cross_vz_r_z / norm_cross_vz_r;
    }

    double dot_vz_r = vz.x * r.x + vz.y * r.y + vz.z * r.z;
    double norm_vz = sqrt(vz.x * vz.x + vz.y * vz.y + vz.z * vz.z);
    double norm_r = sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
    double cosPhi = dot_vz_r / (norm_vz * norm_r + 1e-12);
    double phi = acos(fmin(fmax(cosPhi, -1.0), 1.0));

    double rx = gamma.x;
    double ry = gamma.y;
    double rz = gamma.z;
    double d = 1 - cos(phi);

    newmat3 Q_star;

    Q_star[0][0] = cos(phi) + rx * rx * d;
    Q_star[0][1] = rx * ry * d - rz * sin(phi);
    Q_star[0][2] = rx * rz * d + ry * sin(phi);
    Q_star[1][0] = ry * rx * d + rz * sin(phi);
    Q_star[1][1] = cos(phi) + ry * ry * d;
    Q_star[1][2] = ry * rz * d - rx * sin(phi);
    Q_star[2][0] = rz * rx * d - ry * sin(phi);
    Q_star[2][1] = rz * ry * d + rx * sin(phi);
    Q_star[2][2] = cos(phi) + rz * rz * d;

    double detQ = Q_star[0][0] * (Q_star[1][1] * Q_star[2][2] - Q_star[1][2] * Q_star[2][1]) - Q_star[0][1] * (Q_star[1][0] * Q_star[2][2] - Q_star[1][2] * Q_star[2][0]) + Q_star[0][2] * (Q_star[1][0] * Q_star[2][1] - Q_star[1][1] * Q_star[2][0]);

    if (detQ == 1)
    {
        for (int i = 0; i < DIMENSION; i++)
            for (int j = 0; j < DIMENSION; j++)
                Q[i][j] = Q_star[j][i];
    }
    else
    {
        newmat3 mat = {{1, 0, 0}, {0, 1, 0}, {0, 0, -1}};
        newmat3 QT ;
        transpose(Q_star, QT);
        matrixMultiply(mat, QT, Q);
    }
}


void omega_liutex(const vector u, vector liutex, scalar omega_r)
{

    foreach ()
    {
        double JJ[DIMENSION][DIMENSION];

        scalar s = u.x;
        int i = 0;
        foreach_dimension()
            JJ[0][i++] = center_gradient(s);
        s = u.y; i = 0;
        foreach_dimension()
            JJ[1][i++] = center_gradient(s);
        s = u.z; i = 0;
        foreach_dimension()
            JJ[2][i++] = center_gradient(s);

        /*
                JJ[0][0] = (u.x[1] - u.x[-1]) / (2. * Delta);
                JJ[0][1] = (u.x[0, 1] - u.x[0, -1]) / (2. * Delta);
                JJ[0][2] = (u.x[0, 0, 1] - u.x[0, 0, -1]) / (2. * Delta);
                JJ[1][0] = (u.y[1] - u.y[-1]) / (2. * Delta);
                JJ[1][1] = (u.y[0, 1] - u.y[0, -1]) / (2. * Delta);
                JJ[1][2] = (u.y[0, 0, 1] - u.y[0, 0, -1]) / (2. * Delta);
                JJ[2][0] = (u.z[1] - u.z[-1]) / (2. * Delta);
                JJ[2][1] = (u.z[0, 1] - u.z[0, -1]) / (2. * Delta);
                JJ[2][2] = (u.z[0, 0, 1] - u.z[0, 0, -1]) / (2. * Delta);
        */

        double aa = -(JJ[0][0] + JJ[1][1] + JJ[2][2]);
        // tensor a = gradU[cell];
        // tensor tt = a & a;
        double bb = JJ[0][0] * JJ[1][1] + JJ[0][0] * JJ[2][2] + JJ[1][1] * JJ[2][2] - JJ[0][1] * JJ[1][0] - JJ[1][2] * JJ[2][1] - JJ[2][0] * JJ[0][2];
        double cc = -JJ[0][0] * JJ[1][1] * JJ[2][2] - JJ[0][1] * JJ[1][2] * JJ[2][0] - JJ[0][2] * JJ[2][1] * JJ[1][0] + JJ[0][0] * JJ[1][2] * JJ[2][1] + JJ[1][1] * JJ[2][0] * JJ[0][2] + JJ[2][2] * JJ[0][1] * JJ[1][0];
        double delta = 18.0 * aa * bb * cc - 4.0 * pow(aa, 3) * cc + pow(aa, 2) * pow(bb, 2) - 4.0 * pow(bb, 3) - 27.0 * pow(cc, 2);
        // double qq = pow(aa, 2) * pow(3 * bb, 3) / 9.0;
       
        double qq = (pow(aa, 2) - 3 * bb) / 9.0;
        double rr = (2.0 * pow(aa, 3) - 9.0 * aa * bb + 27.0 * cc) / 54.0;
        delta = -delta / 108.0;
        //double delta = qq * qq * qq  - rr * rr ;

        vec3 Rotex = {0, 0, 0};
        double omega_v=0.0;

            if (delta >0.0)
        {
            double aaaa = -sign(rr) * pow(abs(rr) + sqrt(delta), 1.0 / 3.0);
            double bbbb = 0;

            if (aaaa != 0.0)
                bbbb = qq / aaaa;

             double eig1c_r = -0.5 * (aaaa + bbbb) - aa / 3.0;
            //double eig1c_i = 0.5 * sqrt(3.0) * (aaaa - bbbb);
           // double eig2c_r = -0.5 * (aaaa + bbbb) - aa / 3.0;
            // double eig2c_i = -0.5 * sqrt(3.0) * (aaaa - bbbb);
            double eig3r = aaaa + bbbb - aa / 3.0;
            double delta1 = 0.0, delta2 = 0.0, delta3 = 0.0;
            delta1 = (JJ[0][0] - eig3r) * (JJ[1][1] - eig3r) - JJ[1][0] * JJ[0][1];
            delta2 = (JJ[1][1] - eig3r) * (JJ[2][2] - eig3r) - JJ[1][2] * JJ[2][1];
            delta3 = (JJ[0][0] - eig3r) * (JJ[2][2] - eig3r) - JJ[2][0] * JJ[0][2];

            if (delta1 == 0.0 && delta2 == 0.0 && delta3 == 0.0)
            {
                fprintf(stderr, "delta1-3 are: %f, %f, %f\n", delta1, delta2, delta3);
                exit(1);
            }
            vec3 vr = {0, 0, 0};

            if (fabs(delta1) >= fabs(delta2) && fabs(delta1) >= fabs(delta3))
            {
                vr.x = (-(JJ[1][1] - eig3r) * JJ[0][2] + JJ[0][1] * JJ[1][2]) / delta1;
                vr.y = (JJ[1][0] * JJ[0][2] - (JJ[0][0] - eig3r) * JJ[1][2]) / delta1;
                vr.z = 1.0;
            }
            else if (fabs(delta2) >= fabs(delta1) && fabs(delta2) >= fabs(delta3))
            {
                vr.x = 1.0;
                vr.y = (-(JJ[2][2] - eig3r) * JJ[1][0] + JJ[1][2] * JJ[2][0]) / delta2;
                vr.z = (JJ[2][1] * JJ[1][0] - (JJ[1][1] - eig3r) * JJ[2][0]) / delta2;
            }
            else if (fabs(delta3) >= fabs(delta1) && fabs(delta3) >= fabs(delta2))
            {
                vr.x = (-(JJ[2][2] - eig3r) * JJ[0][1] + JJ[0][2] * JJ[2][1]) / delta3;
                vr.y = 1.0;
                vr.z = (JJ[2][0] * JJ[0][1] - (JJ[0][0] - eig3r) * JJ[2][1]) / delta3;
            }
            else
            {
                fprintf(stderr, "vr error\n");
                exit(1);
            }

            vr.x = vr.x / sqrt(vr.x * vr.x + vr.y * vr.y + vr.z * vr.z);
            vr.y = vr.y / sqrt(vr.x * vr.x + vr.y * vr.y + vr.z * vr.z);
            vr.z = vr.z / sqrt(vr.x * vr.x + vr.y * vr.y + vr.z * vr.z);
            newmat3 Q;
            getOmegaRTransformMatrix(vr, Q);
            // tensor qqq = rotation(z0, vr);
            //  tensor vg = qqq.T() & a;
            // vg = vg & qqq;
            newmat3 localgradU, QJ, QTranspose;
            /*
            for (int i = 0; i < DIMENSION; i++)
                    for (int j = 0; j < DIMENSION; j++)
                    {
                            localgradU[i][j] = 0.0;
                            for (int k = 0; k < DIMENSION; k++)
                                    localgradU[i][j] += Q[i][k] * JJ[k][j];
                    }*/

            matrixMultiply(Q, JJ, QJ);

            // 第二步：计算Q的转置
            transpose(Q, QTranspose);

            // 第三步：计算(QJ)Q^T
            matrixMultiply(QJ, QTranspose, localgradU);

            double alpha = 0.5 * sqrt(pow(localgradU[0][0] - localgradU[1][1], 2) + pow(localgradU[1][0] + localgradU[0][1], 2));
            double beta = 0.5 * (localgradU[1][0] - localgradU[0][1]);
           double b0 = 1e-3;
           double epsilon = b0 * fmax(beta * beta - alpha * alpha, 0.0);
           //omega_v = beta * beta / (beta * beta + alpha * alpha + epsilon);
                  // if ((localgradU[0][0]*vr.x+localgradU[1][0]*vr.y+localgradU[2][0]*vr.z+localgradU[0][1]*vr.x+localgradU[1][1]*vr.y+localgradU[2][1]*vr.z+localgradU[0][2]*vr.x+localgradU[1][2]*vr.y+localgradU[2][2]*vr.z) > 0.0){
            omega_v = beta * beta / (beta * beta + alpha * alpha + eig1c_r * eig1c_r + 0.5 * eig3r * eig3r + epsilon);
           // omega_r[] = 1;

           if (pow(beta, 2) > pow(alpha, 2))
           {
               double rm = 0.0;
               if (beta > 0.0)
                   rm = 2.0 * (beta - alpha);
               else
                   rm = 2.0 * (beta + alpha);

               Rotex.x = rm * vr.x;
               Rotex.y = rm * vr.y;
               Rotex.z = rm * vr.z;
            }
        }


        omega_r[] =  omega_v ;
       

        liutex.x[] = Rotex.x;
        liutex.y[] = Rotex.y;
        liutex.z[] = Rotex.z;
        //if ((Rotex.x*Rotex.x+Rotex.y*Rotex.y+Rotex.z*Rotex.z) > 0.0)
        //if ((localgradU[0][0]*vr.x+localgradU[1][0]*vr.y+localgradU[2][0]*vr.z+localgradU[0][1]*vr.x+localgradU[1][1]*vr.y+localgradU[2][1]*vr.z+localgradU[0][2]*vr.x+localgradU[1][2]*vr.y+localgradU[2][2]*vr.z) > 0.0)
       // omega_r[] =  omega_v ;
        
    }
}
