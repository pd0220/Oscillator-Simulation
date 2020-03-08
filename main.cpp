// oscillator Monte Carlo simulation with given temperature

// including used headers and libraries
#include <iostream>
#include <math.h>
#include <numeric>
#include <random>
#include <string>
#include <sstream>

// lambda to calculate squared values
auto sq = [](auto const a) { return a * a; };

// choose direction (left or right)
auto Direction(std::mt19937 gen)
{
    // uniform distribution in closed interval [0,1]
    std::uniform_real_distribution<> distr(0, 1);
    return distr(gen);
}

// generate real random numbers in closed interval [0,1]
auto UniformRand(std::mt19937 gen)
{
    std::uniform_real_distribution<> distr(0., 1.);
    return distr(gen);
}

// calculate beta * delta energy
// val --> beta * k * a^2
template <typename T>
auto BetaDeltaE(T val, int nPrev, int nNext)
{
    if (nPrev == nNext)
    {
        std::cout << "ERROR: wrong numbers are given" << std::endl;
        std::exit(-1);
    }
    using R = decltype(val * nPrev);
    return (R)0.5 * val * (sq(nNext) - sq(nPrev));
}

// calculate rate
template <typename T>
auto Rate(T val, int nPrev, int nNext)
{
    auto exponent = BetaDeltaE(val, nPrev, nNext);
    // to lower energy level
    if (exponent < 0)
        return 1.;
    // to higher energy level
    else
        return std::exp(-exponent);
}

// main function
// first argument: number of steps in the simulation
// second argument: inital state
// third argument: given beta * k * a^2 (crucial --> defines the system parameters)
int main(int argc, char **argv)
{
    // setup
    if (argc < 4)
    {
        std::cout << "ERROR: not enough argument." << std::endl;
        std::exit(-1);
    }
    // read arguments
    std::string NArg = argv[1], nArg = argv[2], crucialArg = argv[3];
    std::stringstream NStream(NArg), nStream(nArg), crucialStream(crucialArg);
    // define system variables
    int N{0}, nInit{0};
    double crucial{0};
    // first argument
    NStream >> N;
    // second argument
    nStream >> nInit;
    // third argument
    crucialStream >> crucial;

    // random number generation
    std::random_device rd{};
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrInt(0, 1);
    std::uniform_real_distribution<> distrReal(0., 1.);

    // simulation
    // inital state
    int nPrev{nInit};
    for (int t{0}; t < N; t++)
    {
        // next state
        int nNext{0};
        // choose left or rigth direction
        if (distrInt(gen))
            nNext = nPrev + 1;
        else
            nNext = nPrev - 1;

        // calculate rate in chosen direction
        double rate = Rate(crucial, nPrev, nNext);
        if (rate > distrReal(gen))
            nPrev = nNext;
    }

    std::cout << nPrev << std::endl;
}
