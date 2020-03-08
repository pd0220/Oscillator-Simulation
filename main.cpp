// oscillator Monte Carlo simulation with given temperature

// including used headers and libraries
#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include <string>
#include <sstream>
#include <vector>
#include <numeric>

// ------------------------------------------------------------------------------------------------------------------------------------------

// lambda to calculate squared values
auto sq = [](auto const a) { return a * a; };

// ------------------------------------------------------------------------------------------------------------------------------------------

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

// ------------------------------------------------------------------------------------------------------------------------------------------

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

// ------------------------------------------------------------------------------------------------------------------------------------------

// estimation of relaxation time --> first zero value in position
template <typename T1, typename T2>
auto tauEstimate(std::vector<T1> const &time, std::vector<T2> const &position)
{
    // check if we are already in 0
    if (position[0] == 0)
    {
        std::cout << "Simulation started at 0." << std::endl;
        return 0;
    }

    // check for errors
    if (position.size() != time.size())
    {
        std::cout << "ERROR: Container sizes do not match." << std::endl;
        std::exit(-1);
    }

    // determine relaxation time
    int tau{0};
    for (int i{1}; i < static_cast<int>(time.size()); i++)
    {
        if (position[i] == 0)
        {
            tau = time[i];
            return tau;
        }
    }

    // if nothing happens give error message
    std::cout << "ERROR: Relaxation time cannot be estimated.\nConsider rerunning the simulation." << std::endl;
    std::exit(-1);
}

// ------------------------------------------------------------------------------------------------------------------------------------------

// calculate mean of values
template <typename T>
auto Mean(std::vector<T> const &vec)
{
    return std::accumulate(vec.begin(), vec.end(), 0.) / static_cast<double>(vec.size());
}

// ------------------------------------------------------------------------------------------------------------------------------------------

// calculate mean of squared values
template <typename T>
auto MeanSq(std::vector<T> const &vec)
{
    auto addSq = [](auto a, auto b) { return a + sq(b); };
    return std::accumulate(vec.begin(), vec.end(), 0., addSq) / static_cast<double>(vec.size());
}

// ------------------------------------------------------------------------------------------------------------------------------------------

// calculate square of every element in a vector
template <typename T>
auto vecSq(std::vector<T> vec)
{
    int size = static_cast<int>(vec.size());
    std::vector<T> vecSq(size, 0.);
    for (int i{0}; i < size; i++)
    {
        vecSq[i] = sq(vec[i]);
    }
    return vecSq;
}

// ------------------------------------------------------------------------------------------------------------------------------------------

// estimate error
template <typename T1, typename T2>
auto Variance(std::vector<T1> averages, T2 estimator)
{
    using R = decltype(averages[0] + estimator);
    // calculation
    auto add_square = [estimator](R sum, R i) {auto d = i - estimator; return sum + d * d; };
    return std::accumulate(averages.begin(), averages.end(), 0.0, add_square) / static_cast<int>(averages.size());
}

// main function
// first argument: number of steps in the simulation
// second argument: inital state
// third argument: given beta * k * a^2 (crucial --> defines the system parameters)
// fourth argument: file name to save raw data
int main(int argc, char **argv)
{
    // setup
    if (argc < 6)
    {
        std::cout << "ERROR: not enough argument." << std::endl;
        std::exit(-1);
    }

    // read arguments
    std::string NArg = argv[1], nArg = argv[2], crucialArg = argv[3], fileName = argv[4];
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

    // ------------------------------------------------------------------------------------------------------------------------------------------

    // random number generation
    std::random_device rd{};
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrInt(0, 1);
    std::uniform_real_distribution<> distrReal(0., 1.);

    // ------------------------------------------------------------------------------------------------------------------------------------------

    // write to file
    std::ofstream data;
    data.open(fileName);
    data << "time "
         << "state" << std::endl;

    // ------------------------------------------------------------------------------------------------------------------------------------------

    // containers for further analysis
    // store time
    std::vector<int> time(N, 0);
    // store position
    std::vector<int> position(N, 0);

    // ------------------------------------------------------------------------------------------------------------------------------------------

    // simulation
    // inital state
    int nPrev{nInit};
    // initial values to containers
    position[0] = nPrev;
    for (int t{1}; t < N; t++)
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
        // decide to make transition or stay in place
        if (rate > distrReal(gen))
            nPrev = nNext;

        // write step data to file
        data << t << " " << nPrev << std::endl;

        // update conatiners
        time[t] = t, position[t] = nPrev;
    }

    data.close();

    // ------------------------------------------------------------------------------------------------------------------------------------------

    // estimate relaxation time
    int tau = tauEstimate(time, position);

    // delete dynamics from containers during relaxation time
    time.erase(time.begin(), time.begin() + tau);
    position.erase(position.begin(), position.begin() + tau);

    // ------------------------------------------------------------------------------------------------------------------------------------------

    // calulate thermodynamic averages
    double estimator_x = Mean(position), estimator_xSq = MeanSq(position);
    // estimate errors
    double sigma_x = std::sqrt(Variance(position, estimator_x)), sigma_xSq = std::sqrt(Variance(vecSq(position), estimator_xSq));

    std::cout << "E(x): " << estimator_x << " +/- " << sigma_x << std::endl;
    std::cout << "E(x^2): " << estimator_xSq << " +/- " << sigma_xSq << std::endl;
}
