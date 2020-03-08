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

// main function
// first argument: number of steps in the simulation
// second argument: inital state
// third argument: given beta * k * a^2 (crucial --> defines the system parameters)
// fourth argument: file name to save raw data
// fifth argument: number of chunks to analyize
int main(int argc, char **argv)
{
    // setup
    if (argc < 6)
    {
        std::cout << "ERROR: not enough argument." << std::endl;
        std::exit(-1);
    }

    // read arguments
    std::string NArg = argv[1], nArg = argv[2], crucialArg = argv[3], fileName = argv[4], chunksArg = argv[5];
    std::stringstream NStream(NArg), nStream(nArg), crucialStream(crucialArg), chunksStream(chunksArg);
    // define system variables
    int N{0}, nInit{0}, chunks{0};
    double crucial{0};
    // first argument
    NStream >> N;
    // second argument
    nStream >> nInit;
    // third argument
    crucialStream >> crucial;
    // fifth argument
    chunksStream >> chunks;

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

    // new size of containers (in equilibrium)
    int N_eq = static_cast<int>(position.size());

// ------------------------------------------------------------------------------------------------------------------------------------------

    // split to chunks and calculate averages (like jackknife method)
    double numOfValsInChunk = std::floor(N_eq / chunks);

    // containers for thermodynamic averages
    std::vector<double> means(chunks, 0.0);
    std::vector<double> meanSqs(chunks, 0.0);

// ------------------------------------------------------------------------------------------------------------------------------------------

    // calculate averages for chunks
    for (int i{0}; i < chunks - 1; i++)
    {
        // temporary vector to calculate averages
        std::vector<int> tmp(numOfValsInChunk, 0);
        for (int j{0}; j < numOfValsInChunk; j++)
        {
            tmp[j] = position[i * numOfValsInChunk + j];
        }

        // calculate averages
        means[i] = Mean(tmp);
        meanSqs[i] = MeanSq(tmp);

        std::cout << means[i] << " " << meanSqs[i] << std::endl;
    }

    // last chunk might be longer
    double valsInLastChunk = N_eq - (chunks - 1) * numOfValsInChunk;
    std::vector<int> tmp(valsInLastChunk, 0);
    for (int j{0}; j < valsInLastChunk; j++)
    {
        tmp[j] = position[(chunks - 1) * numOfValsInChunk + j];
    }

    means[chunks - 1] = Mean(tmp);
    meanSqs[chunks - 1] = MeanSq(tmp);

    std::cout << means[chunks - 1] << " " << meanSqs[chunks - 1] << std::endl;

// ------------------------------------------------------------------------------------------------------------------------------------------
}
