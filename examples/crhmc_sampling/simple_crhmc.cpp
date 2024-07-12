// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2022-2022 Ioannis Iakovidis

// Contributed and/or modified by Ioannis Iakovidis, as part of Google Summer of
// Code 2022 program.

// Licensed under GNU LGPL.3, see LICENCE file

// References
// Yunbum Kook, Yin Tat Lee, Ruoqi Shen, Santosh S. Vempala. "Sampling with
// Riemannian Hamiltonian
// Monte Carlo in a Constrained Space"


#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "generators/order_polytope_generator.h"
#include "diagnostics/effective_sample_size.hpp"
#include "generators/h_polytopes_generator.h"
#include "volume/sampling_policies.hpp"
#include "ode_solvers/ode_solvers.hpp"
#include "preprocess/crhmc/crhmc_input.h"
#include "preprocess/crhmc/crhmc_problem.h"
#include "sampling/random_point_generators.hpp"
#include "sampling/sampling.hpp"
#include "misc/misc.h"
#include "misc/poset.h"
#include "random.hpp"
#include <vector>
#include "random_walks/random_walks.hpp"
#include "generators/known_polytope_generators.h"
#include "helper_functions.hpp"
#include "volume/rotating.hpp"
#include "convex_bodies/orderpolytope.h"
#include "preprocess/inscribed_ellipsoid_rounding.hpp"


template 
<
typename Polytope,
typename Point,
typename RandomNumberGenerator
>
void sample_cdhr (Polytope &P, RandomNumberGenerator &rng, std::list<Point> &randPoints, unsigned int const&N) {
        Point p = P.ComputeInnerBall().first;
        typedef typename CDHRWalk::template Walk
                <
                        Polytope,
                        RandomNumberGenerator
                > walk;

        typedef RandomPointGenerator <walk> RandomPointGenerator;
        PushBackWalkPolicy push_back_policy;

        RandomPointGenerator::apply(P, p, N, 1, randPoints,
                                    push_back_policy, rng);
}

template 
<
typename Polytope,
typename Point,
typename RandomNumberGenerator
>
void sample_aBW (Polytope &P, RandomNumberGenerator &rng, std::list<Point> &randPoints, unsigned int const&N) {
        
        std::cout << "finding the inner ball" << std::endl;
        Point p = P.ComputeInnerBall().first;
        std::cout << "finished finding the inner ball" << std::endl;
        typedef typename AcceleratedBilliardWalk::template Walk
                <
                        Polytope,
                        RandomNumberGenerator
                > walk;
        typedef RandomPointGenerator <walk> RandomPointGenerator;
        PushBackWalkPolicy push_back_policy;

        RandomPointGenerator::apply(P, p, N, 1, randPoints,
                                    push_back_policy, rng);
}

template 
<
typename MT,
typename VT,
typename NT,
typename Polytope,
typename Point,
typename RandomNumberGenerator
>
void sample_gaBW (Polytope &P, RandomNumberGenerator &rng, std::list<Point> &randPoints, unsigned int const&N) {

        std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
        start = std::chrono::high_resolution_clock::now();


        std::cout << "trying to compute the inner ball" << std::endl;

        Point p = P.ComputeInnerBall().first;
        typedef typename GABW::template Walk
                <
                        Polytope,
                        RandomNumberGenerator
                > walk;
        typedef MultivariateGaussianRandomPointGenerator <walk> RandomPointGenerator;
        PushBackWalkPolicy push_back_policy;

        std::cout << "trying to compute the ellipsoid" << std::endl;

        std::tuple<MT, VT, NT> ellipsoid = compute_inscribed_ellipsoid<MT, EllipsoidType::MAX_ELLIPSOID>
        (P.get_mat(), P.get_vec(), p.getCoefficients(), 500, std::pow(10, -6.0), std::pow(10, -4.0));
        const MT E = get<0>(ellipsoid);

        stop = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> total_time = stop - start;
        std::cout << "Preprocessing Done in " << total_time.count() << '\n';

        start = std::chrono::high_resolution_clock::now();

        walk w(P, p, E, rng);
        for (unsigned int i = 0; i < N; ++i)
        {
            w.template apply(P, p, E, 1, rng);
            push_back_policy.apply(randPoints, p);
        }

        stop = std::chrono::high_resolution_clock::now();

        total_time = stop - start;
        std::cout << "Points sampled in " << total_time.count() << '\n';
        //RandomPointGenerator::apply(P, p, E, N, 1, randPoints,
        //                            push_back_policy, rng);
}


template <typename Polytope, typename MT, typename VT>
double get_hash(Polytope const &P) {
    MT A = P.get_mat();
    double sf = 1;
    for(int i = 0; i < P.dimension(); ++i) {
        double s = 0;
        for(auto j:A.col(i))
        s += j;
        sf = sf * (s + 1.2) / 3.3;
    }
    return sf;
}

double duration, maxPsrf;
unsigned int minEss, nr_samples;

template <int simdLen>
void sample_hpoly(int n_samples = 80000,
                  int n_burns = 20000, int dim = 20, int m = 150, int poly_type = 0, int crhmc_walk = 0) {
                    
    using NT = double;
    using Kernel = Cartesian<NT>;
    using Point = typename Kernel::Point;
    using Func = ZeroScalarFunctor<Point>;
    using Grad = ZeroFunctor<Point>;
    using Hess = ZeroFunctor<Point>;
    using PolytopeType = HPolytope<Point>;
    using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
    using MT = PolytopeType::MT;
    typedef boost::mt19937 PolyRNGType;
    using RNG = BoostRandomNumberGenerator<boost::mt19937, NT>;

    /*
    0 = random skinny
    1 = order poly
    2 = birkoff
    3 = simplex
    4 = skinny cube
    */

    RNG rng(dim);
    PolytopeType HP;
    if(poly_type == 0) {
        HP = skinny_random_hpoly<PolytopeType, NT, PolyRNGType>(dim, m, true, NT(20000), 101);
        std::cout << "Sampling from Random Skinny Polytope of dimension " << dim << " and hash  " << get_hash<PolytopeType, MT, VT>(HP) << std::endl;
    } else if(poly_type == 1) {
        HP = random_orderpoly<PolytopeType, NT>(dim, m, 101);
        std::cout << "Sampling from Order Polytope of dimension " << dim << " and hash  "  << get_hash<PolytopeType, MT, VT>(HP) << std::endl;
    } else if(poly_type == 2) {
        unsigned int aux = ceil(sqrt(dim));
        aux += 1;
        HP = generate_birkhoff<PolytopeType>(aux);
        dim = aux * (aux - 2) + 1;
        std::cout << "Sampling from Birkhoff Polytope of dimension " << aux * (aux - 2) + 1 << " and hash  "  << get_hash<PolytopeType, MT, VT>(HP) << std::endl;
    } else if(poly_type == 3) {
        HP = generate_simplex<PolytopeType>(dim, false);
        std::cout << "Sampling from Simplex Polytope of dimension " << dim << " and hash  "  << get_hash<PolytopeType, MT, VT>(HP) << std::endl;
    } else if(poly_type == 4) {
        HP = generate_skinny_cube<PolytopeType>(dim);
        std::cout << "Sampling from Skinny Cube of dimension " << dim << " and hash  "  << get_hash<PolytopeType, MT, VT>(HP) << std::endl;
    }

    std::ofstream poly_out;
    poly_out.open("cpp_poly.txt");
    MT A = HP.get_mat();
    VT b = HP.get_vec();
    for (unsigned int i = 0; i < A.rows(); i++) {
        for (unsigned int j = 0; j < dim; j++) {
            poly_out << A(i, j) << " ";
        }
        poly_out << " " << b(i) << '\n';
    }
    poly_out.flush();
    poly_out.close();
    if(n_samples == 0)
        return;

    //HP.print();
    Func * f = new Func;
    Grad * g = new Grad;
    std::list<Point> PointList;

    
    std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
    start = std::chrono::high_resolution_clock::now();

    if(crhmc_walk == 1) {
        std::cout << "Using CRHMC walk" << std::endl;
        execute_crhmc< PolytopeType, RNG, std::list<Point>, Grad, Func, Hess, CRHMCWalk, simdLen>(
        HP, rng, PointList, 1, n_samples, n_burns, g, f);
    } else if(crhmc_walk == 0) {
        std::cout << "Using aBW walk" << std::endl;
        sample_aBW(HP, rng, PointList, n_samples);
    } else if(crhmc_walk == 2) {
        std::cout << "Using GaBW walk" << std::endl;
        sample_gaBW<MT, VT, NT>(HP, rng, PointList, n_samples);
    } else if(crhmc_walk == 4) {
        std::pair<Point, NT> InnerBall = HP.ComputeInnerBall();
        std::tuple<MT, VT, NT> res = inscribed_ellipsoid_rounding<MT, VT, NT>(HP, InnerBall.first);
        sample_aBW(HP, rng, PointList, n_samples);
    }

    stop = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> total_time = stop - start;
    std::cout << "Done in " << total_time.count() << '\n';
    duration = total_time.count();
    nr_samples = PointList.size();

    MT samples = MT(dim, PointList.size());

    int i=0;
    for (std::list<Point>::iterator it = PointList.begin(); it != PointList.end(); ++it){
        samples.col(i) = (*it).getCoefficients();
        i++;
    }
    // std::cout << samples << std::endl;
    unsigned int min_ess;
    NT max_psrf;
    effective_sample_size<NT, VT>(samples, min_ess);
    max_psrf = max_interval_psrf<NT,VT,MT>(samples);
    std::cerr << "min_ess: " << min_ess << '\n';
    std::cerr<<"max_psrf: "<< max_psrf <<"\n";
    minEss = min_ess;
    maxPsrf = max_psrf;
    
    std::ofstream samples_stream;
    samples_stream.open("CRHMC_SIMD_" + std::to_string(simdLen) + "_simplex" + "_samples.txt");
    samples_stream << samples.transpose() << std::endl;
    

    
    
    delete f;
    delete g;

}

template<int simdLen>
void run_main(int n_samples = 80000,
              int n_burns = 20000,
              int dimension = 20, int m = 150, int poly_type = 0, int crhmc_walk = 0) {
    sample_hpoly<simdLen>(n_samples, n_burns, dimension, m, poly_type, crhmc_walk);
}

int main(int argc, char *argv[]) {
    
    if (argc != 8) {
        if(argc == 4) // use this just to get matrix output
        {
            // ./simple_crhmc [dim] [facets] [is_order]
            // atoi(argv[1]) = dim
            // atoi(argv[2]) = facets
            // atoi(argv[3]) = 0 for skinny, 1 for order
            run_main<1>(0,0,atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
            exit(1);
        }
        if(argc == 9) // use this to sample many times and output results
        {
            std::ofstream samples_stream;
            samples_stream.open("output" + std::to_string(atoi(argv[1])) + ".txt");
            
            for(int n = 50; n <= 200; n += 50)
            {
                std::cout << "\n\n" << n << '\n';
                int m = n * 10;
                double dur[5][5];
                double psrf[5][5];
                double ess[5][5];
                double nr_s[5][5];
                for(int a = 0; a <= 4; ++a)
                    for(int b = 0; b <= 4; b+=2) {
                        if(atoi(argv[1]) == 1)
                            run_main<1>(60 * n, 10, n, m, a, b);
                        else if(atoi(argv[1]) == 4)
                            run_main<4>(60 * n, 10, n, m, a, b);
                        else if(atoi(argv[1]) == 8)
                            run_main<8>(60 * n, 10, n, m, a, b);
                        else if(atoi(argv[1]) == 16)
                            run_main<16>(60 * n, 10, n, m, a, b);
                        dur[a][b] = duration;
                        psrf[a][b] = maxPsrf;
                        ess[a][b] = minEss;
                        nr_s[a][b] = nr_samples;
                        samples_stream << a << ' ' << b << ' ' << duration << ' ' << maxPsrf << ' ' << minEss << ' ' << nr_samples << std::endl;
                    }
                samples_stream << "\n\n\n";
                /*
                for(int a = 0; a <= 4; ++a) {
                    std::string command = "./get_poly.sh " + std::to_string(n) + " " + std::to_string(n * 10) + " " + std::to_string(a);
                    system(command.c_str());
                    std::ifstream mat_res;
                    mat_res.open("matlab_results.txt");
                    mat_res >> dur[a][3] >> ess[a][3] >> nr_s[a][3] >> psrf[a][3];
                    mat_res.close();
                }*/
                int x;
                samples_stream << "dim = " << n << ",m = " << n*10  << ",random skinny,order,birkoff,simplex,skinny cube\n";
                x = 0;
                samples_stream << "Volesti aBW,runtime," << dur[0][x] << ',' << dur[1][x] << ',' << dur[2][x] << ',' << dur[3][x] << ',' << dur[4][x] << '\n';
                samples_stream << "\"\",samples," << nr_s[0][x] << ',' << nr_s[1][x] << ',' << nr_s[2][x] << ',' << nr_s[3][x] << ',' << nr_s[4][x] << '\n';
                samples_stream << "\"\",psrf," << psrf[0][x] << ',' << psrf[1][x] << ',' << psrf[2][x] << ',' << psrf[3][x] << ',' << psrf[4][x] << '\n';
                samples_stream << "\"\",ess," << ess[0][x] << ',' << ess[1][x] << ',' << ess[2][x] << ',' << ess[3][x] << ',' << ess[4][x] << '\n';
                samples_stream << "\"\",ess / runtime," << ess[0][x] / dur[0][x] << ',' << ess[1][x] / dur[1][x] << ',' << ess[2][x] / dur[2][x] << ',' << ess[3][x] / dur[3][x] << ',' << ess[4][x] / dur[4][x] << '\n';
                samples_stream << "\"\",ess / samples," << ess[0][x] / nr_s[0][x] << ',' << ess[1][x] / nr_s[1][x] << ',' << ess[2][x] / nr_s[2][x] << ',' << ess[3][x] / nr_s[3][x] << ',' << ess[4][x] / nr_s[4][x] << '\n';
                
                samples_stream << "\"\"\n";

                x = 1;
                samples_stream << "Volesti CRHMC,runtime," << dur[0][x] << ',' << dur[1][x] << ',' << dur[2][x] << ',' << dur[3][x] << ',' << dur[4][x] << '\n';
                samples_stream << "\"\",samples," << nr_s[0][x] << ',' << nr_s[1][x] << ',' << nr_s[2][x] << ',' << nr_s[3][x] << ',' << nr_s[4][x] << '\n';
                samples_stream << "\"\",psrf," << psrf[0][x] << ',' << psrf[1][x] << ',' << psrf[2][x] << ',' << psrf[3][x] << ',' << psrf[4][x] << '\n';
                samples_stream << "\"\",ess," << ess[0][x] << ',' << ess[1][x] << ',' << ess[2][x] << ',' << ess[3][x] << ',' << ess[4][x] << '\n';
                samples_stream << "\"\",ess / runtime," << ess[0][x] / dur[0][x] << ',' << ess[1][x] / dur[1][x] << ',' << ess[2][x] / dur[2][x] << ',' << ess[3][x] / dur[3][x] << ',' << ess[4][x] / dur[4][x] << '\n';
                samples_stream << "\"\",ess / samples," << ess[0][x] / nr_s[0][x] << ',' << ess[1][x] / nr_s[1][x] << ',' << ess[2][x] / nr_s[2][x] << ',' << ess[3][x] / nr_s[3][x] << ',' << ess[4][x] / nr_s[4][x] << '\n';
                
                samples_stream << "\"\"\n";

                x = 2;
                samples_stream << "Volesti GaBW,runtime," << dur[0][x] << ',' << dur[1][x] << ',' << dur[2][x] << ',' << dur[3][x] << ',' << dur[4][x] << '\n';
                samples_stream << "\"\",samples," << nr_s[0][x] << ',' << nr_s[1][x] << ',' << nr_s[2][x] << ',' << nr_s[3][x] << ',' << nr_s[4][x] << '\n';
                samples_stream << "\"\",psrf," << psrf[0][x] << ',' << psrf[1][x] << ',' << psrf[2][x] << ',' << psrf[3][x] << ',' << psrf[4][x] << '\n';
                samples_stream << "\"\",ess," << ess[0][x] << ',' << ess[1][x] << ',' << ess[2][x] << ',' << ess[3][x] << ',' << ess[4][x] << '\n';
                samples_stream << "\"\",ess / runtime," << ess[0][x] / dur[0][x] << ',' << ess[1][x] / dur[1][x] << ',' << ess[2][x] / dur[2][x] << ',' << ess[3][x] / dur[3][x] << ',' << ess[4][x] / dur[4][x] << '\n';
                samples_stream << "\"\",ess / samples," << ess[0][x] / nr_s[0][x] << ',' << ess[1][x] / nr_s[1][x] << ',' << ess[2][x] / nr_s[2][x] << ',' << ess[3][x] / nr_s[3][x] << ',' << ess[4][x] / nr_s[4][x] << '\n';
                
                samples_stream << "\"\"\n";
                

                x = 3;
                samples_stream << "Matlab CRHMC,runtime," << dur[0][x] << ',' << dur[1][x] << ',' << dur[2][x] << ',' << dur[3][x] << ',' << dur[4][x] << '\n';
                samples_stream << "\"\",samples," << nr_s[0][x] << ',' << nr_s[1][x] << ',' << nr_s[2][x] << ',' << nr_s[3][x] << ',' << nr_s[4][x] << '\n';
                samples_stream << "\"\",psrf," << psrf[0][x] << ',' << psrf[1][x] << ',' << psrf[2][x] << ',' << psrf[3][x] << ',' << psrf[4][x] << '\n';
                samples_stream << "\"\",ess," << ess[0][x] << ',' << ess[1][x] << ',' << ess[2][x] << ',' << ess[3][x] << ',' << ess[4][x] << '\n';
                samples_stream << "\"\",ess / runtime," << ess[0][x] / dur[0][x] << ',' << ess[1][x] / dur[1][x] << ',' << ess[2][x] / dur[2][x] << ',' << ess[3][x] / dur[3][x] << ',' << ess[4][x] / dur[4][x] << '\n';
                samples_stream << "\"\",ess / samples," << ess[0][x] / nr_s[0][x] << ',' << ess[1][x] / nr_s[1][x] << ',' << ess[2][x] / nr_s[2][x] << ',' << ess[3][x] / nr_s[3][x] << ',' << ess[4][x] / nr_s[4][x] << '\n';
            
                samples_stream << "\"\"\n";
                

                x = 4;
                samples_stream << "Volesti GaBW w rounding," << dur[0][x] << ',' << dur[1][x] << ',' << dur[2][x] << ',' << dur[3][x] << ',' << dur[4][x] << '\n';
                samples_stream << "\"\",samples," << nr_s[0][x] << ',' << nr_s[1][x] << ',' << nr_s[2][x] << ',' << nr_s[3][x] << ',' << nr_s[4][x] << '\n';
                samples_stream << "\"\",psrf," << psrf[0][x] << ',' << psrf[1][x] << ',' << psrf[2][x] << ',' << psrf[3][x] << ',' << psrf[4][x] << '\n';
                samples_stream << "\"\",ess," << ess[0][x] << ',' << ess[1][x] << ',' << ess[2][x] << ',' << ess[3][x] << ',' << ess[4][x] << '\n';
                samples_stream << "\"\",ess / runtime," << ess[0][x] / dur[0][x] << ',' << ess[1][x] / dur[1][x] << ',' << ess[2][x] / dur[2][x] << ',' << ess[3][x] / dur[3][x] << ',' << ess[4][x] / dur[4][x] << '\n';
                
                samples_stream << "\"\"" << std::endl;

                /*
                samples_stream << "\n\n";
                samples_stream << "n = " << n << "             order polytope         skinny polytope\n";
                samples_stream << "CRHMC:  | time:         " << dur[1][1] << "             " << dur[0][1] << '\n';
                samples_stream << "100*n   | psrf:         " << psrf[1][1] << "             " << psrf[0][1] << '\n';
                samples_stream << "samples |  ess:         " << ess[1][1] << "                 " << ess[0][1] << '\n';
                samples_stream << '\n';
                samples_stream << "ABW:    | time:         " << dur[1][0] << "             " << dur[0][0] << '\n';
                samples_stream << "1000*n  | psrf:         " << psrf[1][0] << "             " << psrf[0][0] << '\n';
                samples_stream << "samples |  ess:         " << ess[1][0] << "                 " << ess[0][0] << '\n';
                samples_stream << '\n';
                samples_stream << "Matlab|   time:         " << dur[1][2] << "             " << dur[0][2] << '\n';
                samples_stream << "      |samples:         " << psrf[1][2] << "                " << psrf[0][2] << '\n';
                samples_stream << "      |    ess:         " << ess[1][2] << "             " << ess[0][2] << '\n';

                samples_stream << std::endl;
                */
                

            }
            exit(1);
        }
        std::cerr << "Example Usage: ./simple_crhmc "
                    "[simdLen] [n_samples] [n_burns] [dimension] [facets] [if_order_poly] [if_crhmc_walk]\n";
        std::cerr << "i.e.: ./simple_crhmc 4 1000 500 20 150 1 0\n";
        exit(1);
    }
    // std::cerr << "To plot: python3 ../python_utilities/plot_samples.py <CRHMC_SIMD_4_simplex_samples.txt --save"<<"\n";

    if (atoi(argv[1]) == 1) {
        run_main<1>(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]));
    } else if (atoi(argv[1]) == 4) {
        run_main<4>(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]));
    } else if (atoi(argv[1]) == 8) {
        run_main<8>(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]));
    } else if (atoi(argv[1]) == 16) {
        run_main<16>(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]));
    }
    return 0;
}



/*
            std::ifstream getmat;
            getmat.open("M.txt");
            using NT = double;
            using Kernel = Cartesian<NT>;
            using Point = typename Kernel::Point;
            using Func = ZeroScalarFunctor<Point>;
            using Grad = ZeroFunctor<Point>;
            using Hess = ZeroFunctor<Point>;
            using PolytopeType = HPolytope<Point>;
            using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
            using MT = PolytopeType::MT;
            typedef boost::mt19937 PolyRNGType;
            using RNG = BoostRandomNumberGenerator<boost::mt19937, NT>;

            MT samples = MT(60, 3624);
            for(int i = 0; i < 60; ++i)
            {
                for(int j = 0; j < 3624; ++j)
                {
                    NT a;
                    getmat >> a;
                    samples(i,j) = a;
                }
            }
            unsigned int min_ess = 0;
            effective_sample_size<NT, VT>(samples, min_ess);
            std::cout << "and the ess is....: " << min_ess << '\n';

            exit(1);
*/