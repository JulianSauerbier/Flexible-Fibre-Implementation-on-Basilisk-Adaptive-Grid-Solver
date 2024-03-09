/**
The function below fills the scalar $forcing$ with the regularized source term
of magnitude $force$ located at $x0$. This assumes that $forcing$ is "clean",
i.e. that it is zero everywhere (or that if there is already some non-zero terms
the intention is to include them in the forcing).
*/
#define BGHOSTS 2

#define POS_PBC_X(X) ((u.x.boundary[left] != periodic_bc) ? (X) : (((X - (X0 + L0 / 2)) > L0 / 2.) ? (X)-L0 : (X)))
#define POS_PBC_Y(Y) ((u.x.boundary[top] != periodic_bc) ? (Y) : (((Y - (Y0 + L0 / 2)) > L0 / 2.) ? (Y)-L0 : (Y)))
#define POS_PBC_Z(Z) ((u.x.boundary[front] != periodic_bc) ? (Z) : (((Z - (Z0 + L0 / 2)) > L0 / 2.) ? (Z)-L0 : (Z)))


trace 
void generate_lag_stencils_one_fiber(lagfiber *fiber)
{
        for (int i = 0; i < fiber->nlp; i++)
        {
                fiber->nodes[i].stencil.n = 0;
                double delta = (L0 / (1 << grid->maxdepth));
                for (int ni = -2; ni <= 2; ni++)
                {
                        for (int nj = -2; nj <= 2; nj++)
                        {
                        #if dimension < 3
                        Point point = locate(POS_PBC_X(fiber->nodes[i].post1.x + ni * delta), 
                                             POS_PBC_Y(fiber->nodes[i].post1.y + nj * delta));
                        #else
                        for (int nk = -2; nk <= 2; nk++)
                        {
                        Point point = locate(POS_PBC_X(fiber->nodes[i].post1.x + ni * delta),
                                        POS_PBC_Y(fiber->nodes[i].post1.y + nj * delta),
                                        POS_PBC_Z(fiber->nodes[i].post1.z + nk * delta));
                        #endif
                        // printf("%d %d\n", point.level, grid->maxdepth);
                        // if (point.level >= 0 && point.level != grid->maxdepth)
                        //         fprintf(stderr, "Warning: Lagrangian stencil not fully resolved at node %d.\n", i);
                        cache_append(&(fiber->nodes[i].stencil), point, 0);
                        #if _MPI
                                #if dimension < 3
                                        if (ni == 0 && nj == 0)
                                        {
                                #else
                                        if (ni == 0 && nj == 0 && nk == 0)
                                        {
                                #endif
                                        if (point.level >= 0)
                                                fiber->nodes[i].pid = cell.pid;
                                        else
                                                fiber->nodes[i].pid = -1;
                                        }
                        #endif
                        #if dimension == 3
                        }
                        #endif
                        }
                }
        }
}

trace 
void lag2eul(vector forcing, lagfiber *fiber)
{
    double D_s = fiber->D_s;
    for (int i = 0; i < fiber->nlp; i++)
    {
            foreach_cache(fiber->nodes[i].stencil)
            {
                    if (point.level >= 0)
                    {
                            coord dist;
                            double weight = 0.;
                            dist.x = GENERAL_1DIST(x, fiber->nodes[i].post1.x);
                            dist.y = GENERAL_1DIST(y, fiber->nodes[i].post1.y);
                            if (fabs(dist.x) < 2 * Delta && fabs(dist.y) < 2 * Delta)
                            {
                                weight = (1 + cos(.5 * pi * dist.x / Delta)) * (1 + cos(.5 * pi * dist.y / Delta)) / sq(Delta) * D_s;
                            }
                            #if dimension > 2
                            dist.z = GENERAL_1DIST(z, fiber->nodes[i].post1.z);
                            if (fabs(dist.x) < 2 * Delta && fabs(dist.y) < 2 * Delta && fabs(dist.z) < 2 * Delta)
                            {
                                    weight = (1 + cos(.5 * pi * dist.x / Delta)) * (1 + cos(.5 * pi * dist.y / Delta)) * (1 + cos(.5 * pi * dist.z / Delta)) / (cube(Delta));
                            }
                            #endif
                            foreach_dimension() {
                                    forcing.x[] += DensityD * weight * fiber->nodes[i].LagForce.x;
                                    // fprintf(stderr, "%g %d %g %g %g %g %g %g\n", t, i, fiber->nodes[i].post1.x, dist.x, fiber->nodes[i].LagForce.x, forcing.x[], weight, Delta);
                            }// fprintf(stderr, "%g %d %g %g %g %g %g %g \n", t, i, x, y, fiber->nodes[i].post1.x, fiber->nodes[i].post1.y, , dist.y);
                    }
            }
    }
}

/**
The function below interpolates the eulerian velocities onto the nodes of
the Lagrangian fiber.
*/
trace 
void eul2lag(lagfiber *fiber)
{
        for (int ii = 0; ii < fiber->nlp; ii++)
        {
                foreach_dimension() fiber->nodes[ii].Uib.x = 0.;
                foreach_cache(fiber->nodes[ii].stencil)
                {
                        if (point.level >= 0)
                        {
                                coord dist;
                                double weight = 0.;
                                dist.x = GENERAL_1DIST(x, fiber->nodes[ii].post1.x);
                                dist.y = GENERAL_1DIST(y, fiber->nodes[ii].post1.y);
                                if (fabs(dist.x) < 2 * Delta && fabs(dist.y) < 2 * Delta)
                                {
                                        weight = (1 + cos(.5 * pi * dist.x / Delta)) * (1 + cos(.5 * pi * dist.y / Delta)) / 16.;
                                }
                                #if dimension > 2
                                dist.z = GENERAL_1DIST(z, fiber->nodes[ii].post1.z);
                                if (fabs(dist.x) < 2 * Delta && fabs(dist.y) < 2 * Delta && fabs(dist.z) < 2 * Delta)
                                {
                                        weight = (1 + cos(.5 * pi * dist.x / Delta)) * (1 + cos(.5 * pi * dist.y / Delta)) * (1 + cos(.5 * pi * dist.z / Delta)) / 64.;
                                }
                                #endif
                                foreach_dimension() fiber->nodes[ii].Uib.x += weight * u.x[];
                        }
                }
        }
}

scalar stencils[];
trace 
void tag_ibm_stencils_one_fiber(lagfiber *fiber)
{
        for (int i = 0; i < fiber->nlp; i++)
        {
                foreach()
                {
                        if (point.level >= 0)
                        {
                                coord dist;
                                dist.x = GENERAL_1DIST(x, fiber->nodes[i].post1.x);
                                // Here we don't want the periodic effect on vertical direction
                                dist.y = fabs(y - fiber->nodes[i].post1.y);
                                #if dimension > 2
                                dist.z = GENERAL_1DIST(z, fiber->nodes[i].post1.z);
                                #endif
                                #if dimension < 3
                                if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta)
                                {
                                        stencils[] = sq(dist.x + dist.y) / sq(2. * Delta) * (2. + noise());
                                        // printf("%g %d %g %g\n", t, i, dist.x, dist.y);
                                #else
                                if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta && fabs(dist.z) <= 2 * Delta)
                                {
                                        stencils[] = sq(dist.x + dist.y + dist.z) / cube(2. * Delta) * (2. + noise());
                                #endif
                                }
                                
                        }
                }
        }
}

trace 
void tag_ibm_stencils()
{
        foreach ()
            stencils[] = 0.;
        for (int k = 0; k < NFBS; k++)
                if (fbs.fb[k].isactive)
                        tag_ibm_stencils_one_fiber(&FB(k));
}
