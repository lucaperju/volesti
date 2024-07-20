// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef SAMPLERS_RANDOM_POINT_GENERATORS_HPP
#define SAMPLERS_RANDOM_POINT_GENERATORS_HPP

#include "diagnostics/effective_sample_size.hpp"

template
<
    typename Walk
>
struct RandomPointGenerator
{
    template
    <
        typename Polytope,
        typename Point,
        typename PointList,
        typename WalkPolicy,
        typename RandomNumberGenerator,
        typename Parameters
    >
    static void apply(Polytope& P,
                      Point &p,   // a point to start
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng,
                      Parameters const& parameters)
    {
        Walk walk(P, p, rng, parameters);
        for (unsigned int i=0; i<rnum; ++i)
        {
            walk.template apply(P, p, walk_length, rng);
            policy.apply(randPoints, p);
        }
    }

    template
    <
            typename Polytope,
            typename Point,
            typename PointList,
            typename WalkPolicy,
            typename RandomNumberGenerator
    >
    static void apply(Polytope& P,
                      Point &p,   // a point to start
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng,
                      unsigned int &cnt_iter)
    {
        
        using VT = typename Polytope::VT;
        using MT = typename Polytope::MT;
        bool ok = true; // true if I want to sample until ess = rnum

        std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
        start = std::chrono::high_resolution_clock::now();

        Walk walk(P, p, rng);
        for (unsigned int i=0; ok || i<rnum; ++i)
        {   
            
            stop = std::chrono::high_resolution_clock::now();

            std::chrono::duration<double> total_time = stop - start;
            if(total_time.count() >= 1800)
                break;

            if(ok)
            if(i % 1000 == 0 && i > 0)
            {
                if(i%7 == 0)
                    std::cout << i << std::endl;
                MT samples = MT(P.dimension(), randPoints.size());

                int j=0;
                for (typename PointList::iterator it = randPoints.begin(); it != randPoints.end(); ++it){
                    samples.col(j) = (*it).getCoefficients();
                    j++;
                }
                unsigned int min_ess;
                effective_sample_size<double, VT>(samples, min_ess);
                if(i%7 == 0)
                    std::cout << i << ' ' << min_ess << std::endl;
                if(min_ess >= rnum)
                    break;
            }
            walk.template apply(P, p, walk_length, rng, cnt_iter);
            policy.apply(randPoints, p);
        }
    }
};


template
<
    typename Walk
>
struct MultivariateGaussianRandomPointGenerator
{
    template
    <
        typename Polytope,
        typename Point,
        typename Ellipsoid,
        typename PointList,
        typename WalkPolicy,
        typename RandomNumberGenerator,
        typename Parameters
    >
    static void apply(Polytope& P,
                      Point &p,   // a point to start
                      Ellipsoid const& E,   // ellipsoid representing the Gaussian distribution
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng,
                      Parameters const& parameters)
    {
        Walk walk(P, p, E, rng, parameters);
        for (unsigned int i=0; i<rnum; ++i)
        {
            walk.template apply(P, p, walk_length, rng);
            policy.apply(randPoints, p);
        }
    }

    template
    <
            typename Polytope,
            typename Point,
            typename Ellipsoid,
            typename PointList,
            typename WalkPolicy,
            typename RandomNumberGenerator
    >
    static void apply(Polytope& P,
                      Point &p,   // a point to start
                      Ellipsoid const& E,   // ellipsoid representing the Gaussian distribution
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng,
                      unsigned int &cnt_iter)
    {
        
        using VT = typename Polytope::VT;
        using MT = typename Polytope::MT;
        bool ok = true; // true if I want to sample until ess = rnum

        std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
        start = std::chrono::high_resolution_clock::now();

        Walk walk(P, p, E, rng);
        for (unsigned int i=0; ok || i<rnum; ++i)
        {   
            
            stop = std::chrono::high_resolution_clock::now();

            std::chrono::duration<double> total_time = stop - start;
            if(total_time.count() >= 1800)
                break;

            if(ok)
            if(i % 1000 == 0 && i > 0)
            {
                MT samples = MT(P.dimension(), randPoints.size());

                int j=0;
                for (typename PointList::iterator it = randPoints.begin(); it != randPoints.end(); ++it){
                    samples.col(j) = (*it).getCoefficients();
                    j++;
                }
                unsigned int min_ess;
                effective_sample_size<double, VT>(samples, min_ess);
                if(min_ess >= rnum)
                    break;
            }
            walk.template apply(P, p, walk_length, rng, cnt_iter);
            policy.apply(randPoints, p);
        }
    }
};


template
<
    typename Walk
>
struct GaussianRandomPointGenerator
{
    template
    <
        typename Polytope,
        typename Point,
        typename NT,
        typename PointList,
        typename WalkPolicy,
        typename RandomNumberGenerator
    >
    static void apply(Polytope& P,
                      Point &p,   // a point to start
                      NT const& a_i,
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng)
    {
        Walk walk(P, p, a_i, rng);
        for (unsigned int i=0; i<rnum; ++i)
        {
            walk.template apply(P, p, a_i, walk_length, rng);
            policy.apply(randPoints, p);
        }
    }

    template
    <
            typename Polytope,
            typename Point,
            typename NT,
            typename PointList,
            typename WalkPolicy,
            typename RandomNumberGenerator,
            typename Parameters
    >
    static void apply(Polytope& P,
                      Point &p,   // a point to start
                      NT const& a_i,
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng,
                      Parameters const& parameters)
    {
        Walk walk(P, p, a_i, rng, parameters);

        for (unsigned int i=0; i<rnum; ++i)
        {
            walk.template apply(P, p, a_i, walk_length, rng);
            policy.apply(randPoints, p);
        }
    }
};



template <typename Walk>
struct BoundaryRandomPointGenerator
{
    template
    <
            typename Polytope,
            typename Point,
            typename PointList,
            typename WalkPolicy,
            typename RandomNumberGenerator
    >
    static void apply(Polytope& P,
                      Point &p,   // a point to start
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng)
    {
        Walk walk(P, p, rng);
        Point p1(P.dimension()), p2(P.dimension());
        for (unsigned int i=0; i<rnum; ++i)
        {
            walk.template apply(P, p1, p2, walk_length, rng);
            policy.apply(randPoints, p1);
            policy.apply(randPoints, p2);
        }
    }
};


template
<
    typename Walk
>
struct LogconcaveRandomPointGenerator
{

    template
    <   typename PointList,
        typename WalkPolicy,
        typename RandomNumberGenerator
    >
    static void apply(unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng,
                      Walk &walk)
    {
        typedef double NT;

        for (unsigned int i = 0; i < rnum; ++i)
        {
            // Gather one sample
            walk.apply(rng, walk_length);

            // Use PushBackWalkPolicy
            policy.apply(randPoints, walk.x);
        }
    }
};

template
<
    typename Walk
>
struct CrhmcRandomPointGenerator
{

    template
    <
            typename Polytope,
            typename Point,
            typename PointList,
            typename WalkPolicy,
            typename RandomNumberGenerator,
            typename NegativeGradientFunctor,
            typename NegativeLogprobFunctor,
            typename Parameters
    >
    static void apply(Polytope &P,
                      Point &p,   // a point to start
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng,
                      NegativeGradientFunctor &F,
                      NegativeLogprobFunctor &f,
                      Parameters &parameters,
                      Walk &walk,
                      int simdLen=1,
                      bool raw_output= false)
    {
        typedef typename Walk::MT MT;
        for (unsigned int i = 0; i < std::ceil((float)rnum/simdLen); ++i)
        {
            // Gather one sample
            walk.apply(rng, walk_length);
            if(walk.P.terminate){return;}
            MT x;
            if(raw_output){
              x=walk.x;
            }else{
              x=walk.getPoints();
            }
            if((i + 1) * simdLen > rnum){
              for(int j = 0; j < rnum-simdLen*i; j++){
                Point p = Point(x.col(j));
                policy.apply(randPoints, p);
              }
              break;
            }
            // Use PushBackWalkPolicy
            for(int j=0; j<x.cols();j++){
              Point p = Point(x.col(j));
              policy.apply(randPoints, p);
            }
        }
    }
};

template
<
    typename Walk
>
struct ExponentialRandomPointGenerator
{
    template
    <
        typename Polytope,
        typename Point,
        typename NT,
        typename PointList,
        typename WalkPolicy,
        typename RandomNumberGenerator
    >
    static void apply(Polytope& P,
                      Point &p,   // a point to start
                      Point const& c,   // bias function
                      NT const& T, // temperature/variance
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng)
    {
        Walk walk(P, p, c, T, rng);
        bool success;
        for (unsigned int i=0; i<rnum; ++i)
        {
            success = walk.template apply(P, p, walk_length, rng);
            if (!success) {
                //return;
                throw std::range_error("A generated point is outside polytope");
            }
            policy.apply(randPoints, p);
        }
    }

    template
    <
            typename Polytope,
            typename Point,
            typename NT,
            typename PointList,
            typename WalkPolicy,
            typename RandomNumberGenerator,
            typename Parameters
    >
    static void apply(Polytope& P,
                      Point &p,   // a point to start
                      Point const& c,   // bias function
                      NT const& T, // temperature/variance
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng,
                      Parameters const& parameters)
    {
        Walk walk(P, p, c, T, rng, parameters);
        bool success;

        for (unsigned int i=0; i<rnum; ++i)
        {
            success = walk.template apply(P, p, walk_length, rng);
            if (!success) {
                //return;
                throw std::range_error("A generated point is outside polytope");
            }
            policy.apply(randPoints, p);
        }
    }

};



#endif // SAMPLERS_RANDOM_POINT_GENERATORS_HPP
