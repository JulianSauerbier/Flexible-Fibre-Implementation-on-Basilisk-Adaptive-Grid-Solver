#define ALPHA (-1E5)				//the two large negative constants for the interaction law
#define BETA (-1E2)

bool BC_CLAMP;					//boundary condition at fixed end (true = clamped in condition)
bool GRAVITY;					//gravity enabled or disabled (true = enabled)
coord GCAST = {0., -gCONST};		//gravity vector
double DensityD = 1.; // Rho_fiber / Rho_air density ratio (1 = naturally bouyant case)

typedef struct FNode
{
        coord post0;
        coord post1;
        coord post2;                    // the position of t0, t1, t2 are record as diag position
        coord Xstari;
        double Tt1;                     // Tension force at current time step
        coord Ui;                       // the velocity of the nodes
        coord Uib;                      // the interpolation velocity from fluid to fiber -- later on will be delete
        coord accU;                     // accumulative velocity difference, later on will be delete
        coord LagForce;                 // Lag force
        coord BendForce;                // Bending force
        #if _MPI
        int pid;
        #endif
        Cache stencil;                  // stencil for locating the neighbor points
} FNode;

typedef struct lagfiber
{
        int grid_level;                 // grid level for poisson solver
        int nlp;                       // the number of lag points
        int NN;                        // the number of sections
        double Lr;                      // length of the fiber
        FNode *nodes;                   // the array of lag nodes
        double D_s;                     // Delta_s = the spacing of arclength (assumption)
        bool isactive;

} lagfiber;

typedef struct lagrass
{
        lagfiber fb[NFBS];
        int nfbs;
} lagrass;
lagrass fbs;

#define FB(i) (fbs.fb[i])

//useful macros
#if dimension < 3
#define STENCIL_SIZE 25
#else
#define STENCIL_SIZE 125
#endif

#define ACROSS_PERIODIC(a, b) (fabs(a - b) > L0 / 2.)
#define PERIODIC_1DIST(a, b) (fabs(a - L0 - b) > L0 / 2. ? a + L0 - b : a - L0 - b)
#define GENERAL_1DIST(a, b) (ACROSS_PERIODIC(a, b) ? PERIODIC_1DIST(a, b) : a - b)

#if dimension < 3
#define cnorm(a) (sqrt(sq(a.x) + sq(a.y)))
#define cdot(a, b) (a.x * b.x + a.y * b.y)
#define Vdiff(a, b) ((coord){(a).x - (b).x, (a).y - (b).y})
#define star_Vdiff(a, b) ((coord){(a).x * 2. - (b).x, (a).y * 2. - (b).y})
#define sdot(a) (a.x * a.x + a.y * a.y)
#else
#define cnorm(a) (sqrt(sq(a.x) + sq(a.y) + sq(a.z)))
#define cdot(a, b) (a.x * b.x + a.y * b.y + a.z * b.z)
#define Vdiff(a, b) ((coord){(a).x - (b).x, (a).y - (b).y, (a).z - (b).z})
#define sdot(a) (a.x * a.x + a.y * a.y + a.z * a.z)
#endif

//interaction with fluid (interpolation and coordinate transfer) in:
#include "reg-dirac_flow_v2.h"

//initialises the fibre paramteters
void initialize_fiber(lagfiber *fiber, int nfb)
{       
        double Lr = LFBS;	//length of fiber
        int nlp = NNDS;	//number of nodes
        fiber->nlp = nlp;
        fiber->Lr = Lr;
        fiber->nodes = malloc(nlp * sizeof(FNode));
        fiber->D_s = (double) Lr / (nlp - 1.);
        fiber->isactive = true;
        //initial position
        for (int i = 0; i < nlp; i++)
        {
                fiber->nodes[i].post0.x = X0 + L0/2. + nfb;
                fiber->nodes[i].post0.y = Y0 + L0/2. - Lr + i * fiber->D_s;
                fiber->nodes[i].post1.x = X0 + L0/2. + nfb;
                fiber->nodes[i].post1.y = Y0 + L0/2. - Lr + i * fiber->D_s;
                fiber->nodes[i].Tt1 = 1;

                foreach_dimension(){
                        fbs.fb[nfb].nodes[i].post2.x = 2 * fbs.fb[nfb].nodes[i].post1.x - fbs.fb[nfb].nodes[i].post0.x;
                        fiber->nodes[i].Ui.x = 0.;
                        fiber->nodes[i].accU.x = 0.;
                        fiber->nodes[i].LagForce.x = 0;
                        fiber->nodes[i].BendForce.x = 0;
                }

                fiber->nodes[i].stencil.n = STENCIL_SIZE;
                fiber->nodes[i].stencil.nm = STENCIL_SIZE;
                fiber->nodes[i].stencil.p = malloc(STENCIL_SIZE * sizeof(Index));
        }
         for (int k = 0; k < nlp; k++)
                 fiber->nodes[k].Xstari = star_Vdiff(fbs.fb[nfb].nodes[k].post1, fbs.fb[nfb].nodes[k].post0);  

         while (adapt_wavelet((scalar *){stencils}, (double[]){1.e-5}, maxlevel, minlevel).nf)
        {
                tag_ibm_stencils(&FB(nfb));
                generate_lag_stencils_one_fiber(&FB(nfb));
        }
        BC_CLAMP = false;	//fixed end BC
        GRAVITY = true;		//gravity on/off
}

event init(i = 0)
{
        for (int nfb; nfb < NFBS; nfb++)
                initialize_fiber(&FB(nfb), nfb);
}

void free_fiber(lagfiber *fiber)
{
        for (int i = 0; i < fiber->nlp; i++)
                free(fiber->nodes[i].stencil.p);
        free(fiber->nodes);
}

void free_fibers(lagrass *gfbs)
{
        for (int i = 0; i < gfbs->nfbs; i++)
                if (fbs.fb[i].isactive)
                        free_fiber(&(gfbs->fb[i]));
}

//computes the interaction force
void interactionforce(lagfiber *fiber)
{
        generate_lag_stencils_one_fiber(fiber);
        eul2lag(fiber);
        for(int jj = 0; jj < fiber->nlp-1; jj++)
        {
                foreach_dimension(){
                        fiber->nodes[jj].Ui.x = (fiber->nodes[jj].post1.x - fiber->nodes[jj].post0.x) / dt;
                        fiber->nodes[jj].accU.x += (fiber->nodes[jj].Uib.x - fiber->nodes[jj].Ui.x) * dt;
                        fiber->nodes[jj].LagForce.x = ALPHA * fiber->nodes[jj].accU.x + BETA * (fiber->nodes[jj].Uib.x - fiber->nodes[jj].Ui.x);
                }
        }
}

face vector av[];
vector forcing[];

//acceleration of the fibre based on fluid forces
event acceleration(i++)
{
        foreach ()
                foreach_dimension() forcing.x[] = 0.;
     
        for (int i = 0; i < NFBS; i++)
                if (fbs.fb[i].isactive)
                        lag2eul(forcing, &fbs.fb[i]);
        foreach_face()
                av.x[] = .5 * (forcing.x[] + forcing.x[-1]);
}

//computes bending force of the fibre
void bendingforce(lagfiber *fiber)
{
        int nlp = fiber->nlp;
        double D_s = fiber->D_s;
        double Lgamma = -GAMMA * pow(D_s, -4);

        coord bc_clp = {0., -1.};
        coord Xstari[nlp];
        coord DDXstari[nlp];
        
        for (int k = 0; k < nlp; k++)
                Xstari[k] = fiber->nodes[k].Xstari;

        for (int i = 1; i < nlp-1; i++)
                foreach_dimension()
                        DDXstari[i].x = Xstari[i + 1].x - 2. * Xstari[i].x + Xstari[i - 1].x;

        foreach_dimension()
        {
                fiber->nodes[0].BendForce.x = Lgamma * (DDXstari[2].x - DDXstari[1].x);
                fiber->nodes[1].BendForce.x = Lgamma * (DDXstari[2].x - 2 * DDXstari[1].x);
                if (BC_CLAMP == true)
                        fiber->nodes[nlp - 2].BendForce.x = Lgamma * (2.*bc_clp.x * D_s - 2.*(Xstari[nlp-1].x - Xstari[nlp - 2].x) - 2.*DDXstari[nlp - 2].x + DDXstari[nlp - 3].x);
                else
                        fiber->nodes[nlp - 2].BendForce.x = Lgamma * (-2.*DDXstari[nlp - 2].x + DDXstari[nlp - 3].x); 
        }
        for (int i = 2; i < nlp-2; i++)
                foreach_dimension()
                        fiber->nodes[i].BendForce.x = Lgamma * (DDXstari[i + 1].x - 2. * DDXstari[i].x + DDXstari[i - 1].x);
}

//solves the tension force equation using gauss-seidel method
event Tensionpoisson(i++)
{
        for (int nfb = 0; nfb < NFBS; nfb++)
        {
                dt = dtnext(DT);
                int nlp = fbs.fb[nfb].nlp;
                double D_s = fbs.fb[nfb].D_s;

                interactionforce(&FB(nfb));
                bendingforce(&FB(nfb));

                coord Xstari[nlp];
                coord diffXstari[nlp];
                
                for (int k = 0; k < nlp; k++)
                       Xstari[k] = fbs.fb[nfb].nodes[k].Xstari;
       
                for (int k = 1; k < nlp; k++)
                        diffXstari[k] = Vdiff(Xstari[k], Xstari[k-1]);
                
                // we provide RHS1,2 with i = [0, N-1]; 
                double RHS1[nlp - 1];
                double RHS2[nlp - 1];
                double RHS3[nlp - 1];

                for (int i = 0; i < nlp - 1; i++)
                {
                        double dotXi1_p1 = sdot(Vdiff(fbs.fb[nfb].nodes[i + 1].post1, fbs.fb[nfb].nodes[i].post1));
                        double dotXi1_p0 = sdot(Vdiff(fbs.fb[nfb].nodes[i + 1].post0, fbs.fb[nfb].nodes[i].post0));
                        RHS1[i] = 0.5 / sq(dt) * (1 - 2. * dotXi1_p1 / sq(D_s) + dotXi1_p0 / sq(D_s));
                        RHS2[i] = -1. / sq(D_s) * sdot(Vdiff(fbs.fb[nfb].nodes[i + 1].Ui, fbs.fb[nfb].nodes[i].Ui));
                }
                // RHS3 with i = [1, N-2]
                for (int j = 0; j < nlp - 2; j++)
                {
                        coord diffBendi1 = Vdiff(fbs.fb[nfb].nodes[j + 1].BendForce, fbs.fb[nfb].nodes[j].BendForce);
                        coord diffLagi1  = Vdiff(fbs.fb[nfb].nodes[j + 1].LagForce, fbs.fb[nfb].nodes[j].LagForce);
                        RHS3[j] = -1. / sq(D_s) * cdot(diffXstari[1], Vdiff(diffBendi1, diffLagi1));
                }

                int iteration = 500;
                coord diffBendLagN_1 = Vdiff(fbs.fb[nfb].nodes[nlp - 2].LagForce, fbs.fb[nfb].nodes[nlp - 2].BendForce);
                coord Gra_vector = {GRAVITY * GCAST.x * DensityD, GRAVITY * GCAST.y * DensityD};
                for (int iter = 0; iter < iteration; iter++)
                {
                        // BC at i = 0 and i = N-1
                        fbs.fb[nfb].nodes[0].Tt1 = 1. / 3. / sdot(diffXstari[1]) *
                                                (fbs.fb[nfb].nodes[1].Tt1 * cdot(diffXstari[2], diffXstari[1]) - pow(D_s, 4) * (RHS1[0] + RHS2[0] + RHS3[0]));
                        fbs.fb[nfb].nodes[nlp - 2].Tt1 = 1. / sdot(diffXstari[nlp-1]) *
                                                                (fbs.fb[nfb].nodes[nlp - 3].Tt1 * cdot(diffXstari[nlp-1], diffXstari[nlp-2]) -
                                                                pow(D_s, 4) * (RHS1[nlp - 2] + RHS2[nlp - 2]) + sq(D_s) * cdot(diffXstari[nlp-1], Vdiff(diffBendLagN_1, Gra_vector)));
                        for (int k = 1; k < nlp - 2; k++)
                        {
                                fbs.fb[nfb].nodes[k].Tt1 = 1. / sdot(diffXstari[k]) * (fbs.fb[nfb].nodes[k + 1].Tt1 / 2. * cdot(diffXstari[k + 1], diffXstari[k]) +
                                                        fbs.fb[nfb].nodes[k - 1].Tt1 / 2. * cdot(diffXstari[k], diffXstari[k - 1]) - pow(D_s, 4) / 2. * (RHS1[k] + RHS2[k] + RHS3[k]));
                                
                        }
                }
        }
}

event advectfiber(i++, last)
{       
        for (int nfb = 0; nfb < NFBS; nfb++)
        {
                int nlp = fbs.fb[nfb].nlp;
                double D_s = fbs.fb[nfb].D_s;

                for (int i = 0; i < nlp; i++)
                        foreach_dimension() fbs.fb[nfb].nodes[i].post2.x = 2 * fbs.fb[nfb].nodes[i].post1.x - fbs.fb[nfb].nodes[i].post0.x;
                
                int iteration = 10;
                coord res_X[nlp];
                for (int i = 0; i < nlp-1; i++)
                {
                        foreach_dimension()
                        {
                                res_X[i].x = fbs.fb[nfb].nodes[i].BendForce.x;    // because of the second derivative, the bending force switches sign
                                res_X[i].x += GRAVITY * GCAST.x;
                                res_X[i].x += -fbs.fb[nfb].nodes[i].LagForce.x;
                        }
                }

                for (int iter = 0; iter < iteration; iter++)
                {
                        // When i = 0
                        foreach_dimension() fbs.fb[nfb].nodes[0].post2.x = 1. / (1. / sq(dt) + 2 * fbs.fb[nfb].nodes[0].Tt1 / sq(D_s)) *
                                                                        (2 * fbs.fb[nfb].nodes[0].Tt1 / sq(D_s) * fbs.fb[nfb].nodes[1].post2.x +
                                                                        (2 * fbs.fb[nfb].nodes[0].post1.x - fbs.fb[nfb].nodes[0].post0.x) / sq(dt) +
                                                                        res_X[0].x);
                        for (int j = 1; j < nlp-1; j++)
                        {
                                foreach_dimension()
                                {
                                        fbs.fb[nfb].nodes[j].post2.x = 1. / (sq(D_s / dt) + fbs.fb[nfb].nodes[j].Tt1 + fbs.fb[nfb].nodes[j - 1].Tt1) *
                                                                (fbs.fb[nfb].nodes[j].Tt1 * fbs.fb[nfb].nodes[j + 1].post2.x +
                                                                fbs.fb[nfb].nodes[j - 1].Tt1 * fbs.fb[nfb].nodes[j - 1].post2.x +
                                                                sq(D_s / dt) * (2 * fbs.fb[nfb].nodes[j].post1.x - fbs.fb[nfb].nodes[j].post0.x) +
                                                                sq(D_s) * res_X[j].x);
                                }
                                
                        }
                        
                }
        }
}

//updates the position of every fiber node
event update_fiber(i++, last){
        for (int nfb = 0; nfb < NFBS; nfb++)
        {
                int nlp = fbs.fb[nfb].nlp;
                for (int jj = 0; jj < nlp; jj++)
                {
                        foreach_dimension(){
                                fbs.fb[nfb].nodes[jj].post0.x = fbs.fb[nfb].nodes[jj].post1.x;
                                fbs.fb[nfb].nodes[jj].post1.x = fbs.fb[nfb].nodes[jj].post2.x;
                        }
                }
                for (int k = 0; k < nlp; k++)
                        fbs.fb[nfb].nodes[k].Xstari = star_Vdiff(fbs.fb[nfb].nodes[k].post1, fbs.fb[nfb].nodes[k].post0);
        }
}
