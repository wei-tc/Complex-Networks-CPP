/**
 * Implementation of AdjacencyListGraph.h for CatBrainEdgeList.dat.
 *
 * @author Wei Tao Chi
 * @date 03/09/2018
 */

#include "AdjacencyListGraph.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <numeric>
#include <random>
#include <queue>
#include <set>
#include <map>
#include <ctime>
#include <sstream>

using namespace std;

const int removed = -1;

decomposition::decomposition( unsigned long adjacencyListSize ) : maxCoreness( 0 ),
                                                                  maxLayer( 0 )
{
    cores.resize( adjacencyListSize );
    layers.resize( adjacencyListSize );
}


AdjacencyListGraph::AdjacencyListGraph( const string& fileName )
{
    degreeCount = 0;
    edgeCount = 0;
    nodeCount = 0;
    readAdjacencyList( fileName );
}

AdjacencyListGraph::AdjacencyListGraph( vector<vector<int>> adjacencyList )
        : adjacencyList( std::move( adjacencyList ))
{}

void AdjacencyListGraph::readAdjacencyList( const string& fileName )
{
    ifstream file;
    int nodeA, nodeB;

    setNodeCountAndMappingFromFile( fileName );
    adjacencyList.resize(( unsigned long ) nodeCount );

    // attempt to open file
    file.open( fileName );
    if ( !file.is_open())
    {
        cerr << "failed to open/read file: " << fileName << endl;
        exit( EXIT_FAILURE );
    }

    // skip header lines
    string header;
    for ( int i = 0; i < 4; ++i )
    {
        getline( file, header );
    }

    // read edges
    while ( file >> nodeA >> nodeB )
    {
        int nodeAIdx = mapping.at( nodeA );
        int nodeBIdx = mapping.at( nodeB );

        // -1 since 0-based array, but 1-based node ids
        adjacencyList[nodeAIdx - 1].push_back( nodeBIdx - 1 );

        // uncomment for undirected networks with each edge only appearing once
//        adjacencyList[nodeBIdx - 1].push_back( nodeAIdx - 1 );
        edgeCount++;
    }
    degreeCount = edgeCount * 2;

    file.close();
}

void AdjacencyListGraph::sortNeighboursInAscendingOrderForEachNode
        ( vector<vector<int>>& net )
{
    for ( auto& node : net )
    {
        sort( node.begin(), node.end());
    }
}

int AdjacencyListGraph::setNodeCountAndMappingFromFile( const string& fileName )
{
    ifstream file;
    string line;

    file.open( fileName );
    if ( !file.is_open())
    {
        cerr << "failed to open/read file: " << fileName << endl;
        exit( EXIT_FAILURE );
    }

    // skip header lines
    string header;
    for ( int i = 0; i < 4; ++i )
    {
        getline( file, header );
    }

    int max = 0;
    int node1, node2;
    int nodeIdx = 0;
    while ( file >> node1 >> node2 )
    {
        if ( mapping.find( node1 ) == mapping.end()) // node1 not yet added
        {
            mapping.insert( {node1, ++nodeIdx} );
        }

        if ( mapping.find( node2 ) == mapping.end()) // node2 not yet added
        {
            mapping.insert( {node2, ++nodeIdx} );
        }

        if ( node1 > max )
        {
            max = node1;
        }

        if ( node2 > max )
        {
            max = node2;
        }
    }
    file.close();

    nodeCount = static_cast<int>(mapping.size());

    return max;
}

void AdjacencyListGraph::convertToGnmNetwork()
{
    vector<vector<int>> gnmNetwork(( unsigned long ) nodeCount );

    // used to generate random numbers [0, nodeCount)
    random_device rd;
    mt19937 rng( rd());
    uniform_int_distribution<int> dist( 0, nodeCount - 1 );

    for ( int i = 0; i < edgeCount; ++i )
    {
        int node = dist( rng );
        int neighbour = dist( rng );
        gnmNetwork[node].push_back( neighbour );
        gnmNetwork[neighbour].push_back( node );
    }

    // sort so that duplicates can be easily skipped when calculating the
    // clustering coefficient
    sortNeighboursInAscendingOrderForEachNode( gnmNetwork );
    adjacencyList = gnmNetwork;
}

void AdjacencyListGraph::convertToConfigurationModel()
{
    vector<int> stubs;
    setStubs( stubs );
    shuffleStubs( stubs );
    vector<vector<int>> shuffledNetwork(( unsigned long ) nodeCount );

    for ( int i = 0; i < stubs.size() - 1; i += 2 )
    {
        int node = stubs[i];
        int neighbour = stubs[i + 1];
        shuffledNetwork[node].push_back( neighbour );
        shuffledNetwork[neighbour].push_back( node );
    }

    // sort so that duplicates can be easily skipped when calculating the
    // clustering coefficient
    sortNeighboursInAscendingOrderForEachNode( shuffledNetwork );
    adjacencyList = shuffledNetwork;
}

void AdjacencyListGraph::setStubs( vector<int>& stubs )
{
    // node id start from 1, but set to start from 0,
    // since vectors are 0-indexed
    int i = -1;
    for ( const auto& node : adjacencyList )
    {
        i++;
        for ( int neighbour : node )
        {
            stubs.push_back( i );
        }
    }
}

void AdjacencyListGraph::shuffleStubs( vector<int>& stubs )
{
    shuffle( stubs.begin(), stubs.end(), mt19937( random_device()()));
}

void AdjacencyListGraph::printDegreeDistribution()
{
    vector<int> degreeDistribution = getDegreeDistribution();
    int i = 0;
    for ( double d : degreeDistribution )
    {
        cout << i++ << " degrees" << ": " << d << " node(s)" << endl;
    }
}

vector<int> AdjacencyListGraph::getDegreeDistribution() const
{
    // +1 for where the node has 0 degrees
    vector<int> degreeDistribution( adjacencyList.size(), 0 );

    // get degree frequencies
    for ( auto const& node : adjacencyList )
    {
        degreeDistribution[node.size()]++;
    }

    return degreeDistribution;
}

vector<double> AdjacencyListGraph::getNormalisedDegreeDistribution()
{
    vector<int> degreeDistribution = getDegreeDistribution();
    return normaliseDegreeDistribution( degreeDistribution );
}

vector<double> AdjacencyListGraph::normaliseDegreeDistribution(
        vector<int>& distribution )
{
    vector<double> normalisedDistribution;
    for ( int& numberOfDegreeKNodes : distribution )
    {
        normalisedDistribution.push_back( numberOfDegreeKNodes
                                          / ( double ) nodeCount );
    }

    return normalisedDistribution;
}

double AdjacencyListGraph::getVarianceOfDegreeDistribution()
{
    vector<int> degreeDistribution = getDegreeDistribution();

    auto nodeCount = ( int ) adjacencyList.size();
    double mean = calculateDegreeDistributionMean( degreeDistribution );

    double variance = 0;
    for ( auto& degree: degreeDistribution )
    {
        variance += (pow( degree - mean, 2 ));
    }

    return variance / ( double ) nodeCount;
}

double AdjacencyListGraph::calculateDegreeDistributionMean(
        const std::vector<int>& distribution )
{

    double total = accumulate( distribution.begin(),
                               distribution.end(),
                               0.0 );
    double mean = total / ( double ) nodeCount;
    return mean;
}

string AdjacencyListGraph::getClusteringCoefficient( const AdjacencyListGraph& network )
{
    vector<double> clusteringCoefficients = network.getLocalClusteringCoefficients();
    stringstream ss;

    for ( double clustering : clusteringCoefficients )
    {
        ss << clustering << ",";
    }

    string output = ss.str();
    output.resize( output.length() - 1 ); // to remove extra comma

    return output;
}

void AdjacencyListGraph::printClusteringCoefficients()
{
    vector<double> clusteringCoefficients = getLocalClusteringCoefficients();
    int i = 1; // since node ids start from 1

    for ( double c : clusteringCoefficients )
    {
        cout << "Node " << i++ << ": " << c << endl;
    }
}

vector<double> AdjacencyListGraph::getLocalClusteringCoefficients() const
{
    unsigned long size = adjacencyList.size();

    // max degree is the number of nodes in an undirected graph
    // with max one edge between two nodes
    vector<double> clusteringCoefficients( size, 0.0 );
    int triads;
    int triangles;

    for ( int node = 0; node < size; ++node )
    {
        triads = 0;
        triangles = 0;

        // loop over neighbours
        for ( int nb = 0; nb < adjacencyList[node].size(); ++nb )
        {
            int neighbour = adjacencyList[node][nb];

            // skip reflexive edges (self-loops), which cannot form triads
            bool isReflexiveEdge = (node == neighbour);
            if ( isReflexiveEdge )
            {
                continue;
            }

            // Skip duplicate edges. This assumes that every vector
            // within the adjacency list is sorted.
            bool isDuplicateEdge = (neighbour == adjacencyList[node][nb - 1]);
            if ( nb != 0 && isDuplicateEdge )
            {
                // since nb = 0 is the first neighbour in the adjacency list for
                // that particular node, it cannot be a duplicate edge
                continue;
            }

            // loop over higher neighbours
            for ( int high = nb + 1; high < adjacencyList[node].size(); ++high )
            {
                int higherNeighbour = adjacencyList[node][high];
                triads++;

                if ( find( adjacencyList[neighbour].begin(),
                           adjacencyList[neighbour].end(),
                           higherNeighbour ) != adjacencyList[neighbour].end())
                {
                    triangles++;
                }
            }
        }

        clusteringCoefficients[node] = (triangles == 0 || triads == 0)
                                       ? 0
                                       : triangles / ( double ) triads;
    }

    return clusteringCoefficients;
}

double AdjacencyListGraph::getAverageClusteringCoefficient()
{
    auto localClusteringCoefficients = getLocalClusteringCoefficients();
    double total = accumulate( localClusteringCoefficients.begin(),
                               localClusteringCoefficients.end(),
                               0.0 );
    double average = total / nodeCount;
    return average;
}

void AdjacencyListGraph::printDegreeForEachNodeForExportToR()
{
    cout << "DEGREE OF EACH NODE" << endl;
    for ( const auto& n : adjacencyList )
    {
        cout << n.size() << ", ";
    }
    cout << endl;
}


void AdjacencyListGraph::printDegreeAgainstLocalClusteringCoefficientForEachNode()
{
    vector<double> localClusteringCoefficients = getLocalClusteringCoefficients();

    for ( int i = 0; i < nodeCount; ++i )
    {
        // +1 because nodes are labelled from 1
        cout << i + 1 << ": " << adjacencyList[i].size() << " "
             << localClusteringCoefficients[i] << endl;
    }
}

double AdjacencyListGraph::getAssortativityCoefficient() const
{
    double mean = getAverageExcessDegree();
    double variance = calculateExcessDegreeVariance( mean );

    double excessKiKj = 0;
    for ( int node = 0; node < nodeCount; ++node )
    {
        auto degreeKi = ( int ) adjacencyList[node].size();

        for ( int nb = 0; nb < adjacencyList[node].size(); ++nb )
        {
            int neighbour = adjacencyList[node][nb];
            auto degreeKj = ( int ) adjacencyList[neighbour].size();

            int excessKi = degreeKi - 1;
            int excessKj = degreeKj - 1;
            excessKiKj += (excessKi * excessKj);
        }
    }

    return ((excessKiKj / degreeCount) - pow( mean, 2 )) / variance;
}

double AdjacencyListGraph::getAverageExcessDegree() const
{
    double mean = 0;
    for ( const auto& node : adjacencyList )
    {
        auto degree = ( int ) node.size();
        int excessDegree = degree - 1;
        mean += (degree * excessDegree);
    }

    return mean / ( double ) degreeCount;
}

double AdjacencyListGraph::calculateExcessDegreeVariance( double mean ) const
{
    double variance = 0.0;

    for ( const auto& node : adjacencyList )
    {
        auto degree = static_cast<int>(node.size());
        int excessDegree = degree - 1;
        variance += (pow( excessDegree - mean, 2 ) * degree);
    }

    return variance / ( double ) degreeCount;
}

vector<vector<int>> AdjacencyListGraph::getShortestPathLengthForAllNodes() const
{
    vector<vector<int>> shortestPathLengths;
    for ( int i = 0; i < adjacencyList.size(); ++i )
    {
        // bfs starting from each node
        shortestPathLengths.push_back( getShortestPathLength( i ));
    }

    return shortestPathLengths;
}

vector<int> AdjacencyListGraph::getShortestPathLength( int sourceNode ) const
{
    // BFS
    // step 0a
    vector<int> distances(( unsigned long ) nodeCount );
    fill( distances.begin(), distances.end(), -1 );
    distances[sourceNode] = 0;

    // step 0b
    queue<int> q;
    q.push( sourceNode );

    // step 1
    while ( !q.empty())
    {
        int node = q.front();
        q.pop();

        // step 2
        int distance = distances[node];

        // step 3
        for ( auto neighbour : adjacencyList[node] )
        {
            if ( distances[neighbour] == -1 )
            {
                distances[neighbour] = distance + 1;
                q.push( neighbour );
            }
        }
    }

    return distances;
}

const double AdjacencyListGraph::getAverageShortestPathLength(
        vector<vector<int>> distances ) const
{
    int sum = 0;
    int totalPathCount = 0;

    for ( auto distancesFromStartingNode : distances )
    {
        for ( auto distance : distancesFromStartingNode )
        {
            if ( distance > 0 ) // exclude starting node and unvisited nodes
            {
                sum += distance;
                totalPathCount++;
            }
        }
    }

    return sum / ( double ) totalPathCount;
}

const int AdjacencyListGraph::getDiameterOfNetwork(
        vector<vector<int>> distances ) const
{
    int diameter = -1;

    for ( const auto& distancesFromStartingNode : distances )
    {
        for ( const auto& distance : distancesFromStartingNode )
        {
            if ( distance > diameter )
            {
                diameter = distance;
            }
        }
    }

    return diameter;
}

const vector<double> AdjacencyListGraph::getClosenessCentralityForAllNodes(
        vector<vector<int>> distances ) const
{
    vector<double> closenessCentralities;

    for ( const auto& distancesFromNode : distances )
    {
        closenessCentralities.push_back(
                getClosenessCentrality( distancesFromNode ));
    }

    return closenessCentralities;
}

double AdjacencyListGraph::getClosenessCentrality( vector<int> distances ) const
{
    int sum = 0;
    int numDistances = 0;

    for ( auto distance : distances )
    {
        if ( distance > 0 ) // skip starting node and nodes not in the component
        {
            sum += distance;
            numDistances++;
        }
    }

    return numDistances / ( double ) sum;
}

vector<int> AdjacencyListGraph::selfAvoidingDegreeSeekingLocalSearch( int source,
                                                                      int target ) const
{
    int current = source;
    vector<Colour> colours( static_cast<unsigned long>(nodeCount));
    fill( colours.begin(), colours.end(), white );
    vector<int> localSearch; // includes source and target node

    while ( current != target )
    {
        localSearch.push_back( current );
        setNodeColour( colours, current );

        // visit next node
        current = getTargetOrRandomHighestDegreeWhiteOrGreyNeighbour( current,
                                                                      target,
                                                                      colours );
    }

    localSearch.push_back( current ); // add target node

    return localSearch;
}

vector<int> AdjacencyListGraph::getNeighboursShuffledThenSortedByDegreeByDesc( int current ) const
{
    vector<int> neighboursSorted;

    // get neighbours
    for ( const auto& neighbours : adjacencyList[current] )
    {
        neighboursSorted.push_back( neighbours );
    }

    // shuffle vector so that when this function is called in getRandomHighestDegreeNeighbour
    // and, if there are multiple nodes of the same degree and colour,
    // a random one will be selected as the potential next node
    shuffle( neighboursSorted.begin(),
             neighboursSorted.end(),
             mt19937( random_device()()));

    // stable_sort to ensure randomness is maintained
    stable_sort( neighboursSorted.begin(),
                 neighboursSorted.end(),
                 [this]( int a, int b ) {
                     return adjacencyList[a].size() > adjacencyList[b].size();
                 }
               );

    return neighboursSorted;
}

int AdjacencyListGraph::getTargetOrRandomHighestDegreeWhiteOrGreyNeighbour( int current,
                                                                            int target,
                                                                            vector<Colour>& colours ) const
{
    vector<int> neighboursSorted = getNeighboursShuffledThenSortedByDegreeByDesc( current );
    int greyNeighbour = -1;
    int whiteNeighbour = -1;
    bool isGreyNeighbourSet = false;
    bool isWhiteNeighbourSet = false;
    bool isNeighbourOfTarget = false;

    // get the potential next nodes
    for ( int neighbour : neighboursSorted )
    {
        Colour neighbourColour = colours[neighbour];

        if ( neighbour == target )
        {
            isNeighbourOfTarget = true;
            break;
        }

        if ( !isWhiteNeighbourSet && neighbourColour == white )
        {
            whiteNeighbour = neighbour;
            isWhiteNeighbourSet = true;
            continue;
        }
        if ( !isGreyNeighbourSet && neighbourColour == grey )
        {
            greyNeighbour = neighbour;
            isGreyNeighbourSet = true;
            continue;
        }
    }

    // return the 'best' neighbour of the current node
    if ( isNeighbourOfTarget )
    {
        return target;
    }
    else if ( isWhiteNeighbourSet )
    {
        return whiteNeighbour;
    }
    else
    {
        return greyNeighbour;
    }
}

void AdjacencyListGraph::setNodeColour( std::vector<Colour>& colours,
                                        int current ) const
{
    bool isAllVisited = true;
    for ( const auto& neighbour : adjacencyList[current] )
    {
        if ( colours[neighbour] == white )
        {
            isAllVisited = false;
            break;
        }
    }

    if ( isAllVisited )
    {
        colours[current] = black;
    }
    else
    {
        colours[current] = grey;
    }
}

vector<int> AdjacencyListGraph::degreeBiasedRandomWalk( int source,
                                                        int target )
{
    int current = source;
    vector<int> randomWalk; // includes the source and target node

    while ( current != target )
    {
        randomWalk.push_back( current );

        // visit next node
        bool isNeighbourOfTarget = isCurrentNodeNeighbourOfTarget( current, target );
        current = (isNeighbourOfTarget) ? target
                                        : getRandomHighestDegreeNeighbour( current );
    }

    randomWalk.push_back( current ); // add target node

    return randomWalk;
}

std::vector<unsigned long> AdjacencyListGraph::degreeBiasedRandomWalkByDuration( int source,
                                                                                 time_t duration )
{
    int current = source;
    vector<unsigned long> randomWalk; // excludes source node

    time_t startTime = time( nullptr );
    time_t elapsedTime = 0;
    while ( elapsedTime < duration )
    {
        // visit next random neighbour of highest degree
        current = getRandomHighestDegreeNeighbour( current );
        randomWalk.push_back( static_cast<unsigned long>(current));

        time_t currentTime = time( nullptr );
        elapsedTime = currentTime - startTime;
    }

    return randomWalk;
}

bool AdjacencyListGraph::isCurrentNodeNeighbourOfTarget( int current, int target ) const
{
    for ( const auto& neighbour : adjacencyList[current] )
    {
        if ( neighbour == target )
        {
            return true;
        }
    }

    return false;
}

int AdjacencyListGraph::getRandomHighestDegreeNeighbour( int current )
{
    random_device rd;
    mt19937 rng( rd());

    // if the currentNode's neighbour (idx = 3) has degree 5,
    // neighboursByDegree contains the idx 3, five times
    vector<int> neighboursByDegree = getNeighboursByDegree( current );
    auto neighboursByDegreeSize = static_cast<int>(neighboursByDegree.size());
    uniform_int_distribution<int> dist( 0, neighboursByDegreeSize - 1 );
    int randomNeighbour = dist( rng );

    return neighboursByDegree[randomNeighbour];
}

vector<int> AdjacencyListGraph::getNeighboursByDegree( int current ) const
{
    vector<int> neighboursByDegree;

    for ( const auto& neighbour : adjacencyList[current] )
    {
        auto degree = static_cast<int>(adjacencyList[neighbour].size());
        for ( int i = 0; i < degree; ++i )
        {
            neighboursByDegree.push_back( neighbour );
        }
    }

    return neighboursByDegree;
}

decomposition AdjacencyListGraph::onionAndKCoreDecomposition() const
{
    vector<int> nodes = getNodes();
    vector<int> degreeOfNodes = getDegreeOfEachNode();
    struct decomposition decomp = decomposition( adjacencyList.size());
    int core = 1;
    int layer = 1;
    int nodesRemaining = nodeCount;

    while ( nodesRemaining != 0 )
    {
        vector<int> thisLayer = getLayer( nodes, degreeOfNodes, core );
        for ( int nodeInLayer : thisLayer )
        {
            decomp.cores[nodeInLayer] = core;
            decomp.layers[nodeInLayer] = layer;

            for ( int neighbour : adjacencyList[nodeInLayer] )
            {
                if ( nodes[neighbour] != removed )
                {
                    degreeOfNodes[neighbour]--;
                }
            }

            nodes[nodeInLayer] = removed;
            degreeOfNodes[nodeInLayer] = removed;
            nodesRemaining--;
        }

        layer++;
        int minDegree = getMinDegree( degreeOfNodes );
        if ( minDegree >= core + 1 )
        {
            core = minDegree;
        }
    }

    decomp.maxCoreness = core;
    decomp.maxLayer = layer - 1; // -1 since the last layer has no nodes in the while loop

    return decomp;
}

vector<int> AdjacencyListGraph::getOnionDistribution( vector<int> layers, int maxLayer ) const
{
    vector<int> onionDistribution( static_cast<unsigned long>(maxLayer));

    for ( const auto& layer : layers )
    {
        // -1 since the vector is 0-based, but the layers are 1-based
        onionDistribution[layer - 1]++;
    }

    return onionDistribution;
}

void AdjacencyListGraph::printOnionDistributionForExportToR( const vector<int>& onionDistribution,
                                                             const string& datasetFilename )
{
    cout << "ONION DISTRIBUTION" << endl;
    for ( int layer : onionDistribution )
    {
        cout << layer << ", ";
    }
    cout << endl;
}

string AdjacencyListGraph::getLayerOfEachNodeForExportToR( const vector<int>& layers,
                                                           const string& datasetFilename )
{
    stringstream ss;

    for ( int layer : layers )
    {
        ss << layer << ",";
    }

    string output = ss.str();
    output.resize( output.length() - 1 ); // remove extra comma

    return output;
}

vector<int> AdjacencyListGraph::getKCoreDistribution( vector<int> cores, int maxCoreness ) const
{
    vector<int> kCoreDistribution( static_cast<unsigned long>(maxCoreness));

    for ( const auto& core : cores )
    {
        // -1 since the vector is 0-based, but the k-core is 1-based
        kCoreDistribution[core - 1]++;
    }

    return kCoreDistribution;
}

void AdjacencyListGraph::printCorenessAndLayerOfEachNode( const decomposition& kCD )
{
    cout << "node\tlayer\tcoreness" << endl;
    for ( int i = 0; i < kCD.cores.size(); ++i )
    {
        // i + 1 since the nodes are 1-based
        cout << i + 1 << "\t\t" << kCD.layers[i] << "\t\t" << kCD.cores[i] << endl;
    }
    cout << endl;
}

const string& AdjacencyListGraph::printKCoreDistributionForExportToR(
        const vector<int>& kCoreDistribution, const string& datasetFilename )
{
    cout << "K-CORE DISTRIBUTION" << endl;
    for ( int core : kCoreDistribution )
    {
        cout << core << ", ";
    }
    cout << endl;
}

string AdjacencyListGraph::getCorenessOfEachNodeForExportToR( const vector<int>& cores,
                                                              const string& datasetFilename )
{
    stringstream ss;

    for ( int core : cores )
    {
        ss << core << ",";
    }

    string output = ss.str();
    output.resize( output.length() - 1 ); // remove extra comma

    return output;
}

vector<int> AdjacencyListGraph::getNodes() const
{
    vector<int> nodes( static_cast<unsigned long>(nodeCount));
    iota( nodes.begin(), nodes.end(), 0 );

    return nodes;
}

vector<int> AdjacencyListGraph::getDegreeOfEachNode() const
{
    vector<int> degreeOfNodes;

    for ( const auto& node : adjacencyList )
    {
        auto degree = static_cast<int>(node.size());
        degreeOfNodes.push_back( degree );
    }

    return degreeOfNodes;
}


vector<int> AdjacencyListGraph::getLayer( const vector<int>& nodes,
                                          const vector<int>& degreeOfNeighbours,
                                          int core ) const
{
    vector<int> layer;

    for ( const auto& node : nodes )
    {
        if ( node == removed )
        {
            continue;
        }

        int degree = degreeOfNeighbours[node];
        if ( degree <= core )
        {
            layer.push_back( node );
        }
    }

    return layer;
}

int AdjacencyListGraph::getMinDegree( vector<int> degreeOfNodes ) const
{
    int minDegree = numeric_limits<int>::max();
    for ( int degree : degreeOfNodes )
    {
        if ( degree == removed )
        {
            continue;
        }

        if ( degree < minDegree )
        {
            minDegree = degree;
        }
    }

    // all nodes have been removed, and their coreness and layer determined
    if ( minDegree == numeric_limits<int>::max())
    {
        minDegree = removed;
    }

    return minDegree;
}

const int AdjacencyListGraph::getDegreeCount() const
{
    return degreeCount;
}

const int AdjacencyListGraph::getEdgeCount() const
{
    return edgeCount;
}

const int AdjacencyListGraph::getNodeCount() const
{
    return nodeCount;
}

const vector<vector<int>>& AdjacencyListGraph::getAdjacencyList() const
{
    return adjacencyList;
}

const map<int, int>& AdjacencyListGraph::getMapping() const
{
    return mapping;
}
