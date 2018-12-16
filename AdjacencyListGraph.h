/**
 * Reads in a graph comprised of an adjacency list and computes the degree
 * and local clustering coefficient distributions. Able to convert to G(N,M)
 * and Configuration Model random networks. Able to calculate
 * closeness centrality, average shortest path lengths, and diameter of the
 * network.
 *
 * @author Wei Tao Chi
 * @date 03/09/2018
 */

#ifndef ANALYSEADJACENCYLIST_H
#define ANALYSEADJACENCYLIST_H

#include <vector>
#include <string>
#include <unordered_map>
#include <map>

enum Colour
{
    white, grey, black
};

struct decomposition
{
    std::vector<int> cores;
    int maxCoreness;
    std::vector<int> layers;
    int maxLayer;

    explicit decomposition( unsigned long adjacencyListSize );
};


class AdjacencyListGraph
{
public:
    explicit AdjacencyListGraph( const std::string& fileName );

    explicit AdjacencyListGraph( std::vector<std::vector<int>> adjacencyList );

    int setNodeCountAndMappingFromFile( const std::string& fileName );

    void convertToGnmNetwork();

    void convertToConfigurationModel();

    void printDegreeDistribution();

    std::vector<int> getDegreeDistribution() const;

    std::vector<double> getNormalisedDegreeDistribution();

    double getVarianceOfDegreeDistribution();

    void printClusteringCoefficients();

    std::string getClusteringCoefficient( const AdjacencyListGraph& network );

    std::vector<double> getLocalClusteringCoefficients() const;

    double getAverageClusteringCoefficient();

    void printDegreeForEachNodeForExportToR();

    void printDegreeAgainstLocalClusteringCoefficientForEachNode();

    double getAssortativityCoefficient() const;

    std::vector<std::vector<int>> getShortestPathLengthForAllNodes()
    const;

    std::vector<int> getShortestPathLength( int sourceNode ) const;

    const double getAverageShortestPathLength(
            std::vector<std::vector<int>> distances ) const;

    const int getDiameterOfNetwork(
            std::vector<std::vector<int>> distances ) const;

    const std::vector<double> getClosenessCentralityForAllNodes
            ( std::vector<std::vector<int>> distances ) const;

    std::vector<int> selfAvoidingDegreeSeekingLocalSearch( int source,
                                                           int target ) const;

    std::vector<int> degreeBiasedRandomWalk( int source, int target );

    std::vector<unsigned long>
    degreeBiasedRandomWalkByDuration( int source, time_t duration );

    decomposition onionAndKCoreDecomposition() const;

    std::vector<int> getKCoreDistribution( std::vector<int> cores, int maxCoreness ) const;

    std::vector<int> getOnionDistribution( std::vector<int> layers, int maxLayer ) const;

    void printOnionDistributionForExportToR( const std::vector<int>& onionDistribution,
                                             const std::string& datasetFilename );

    std::string getLayerOfEachNodeForExportToR( const std::vector<int>& layers,
                                                const std::string& datasetFilename );

    void printCorenessAndLayerOfEachNode( const decomposition& kCD );

    const std::string&
    printKCoreDistributionForExportToR( const std::vector<int>& kCoreDistribution,
                                        const std::string& datasetFilename );

    std::string getCorenessOfEachNodeForExportToR( const std::vector<int>& cores,
                                                   const std::string& datasetFilename );

    std::vector<int> getDegreeOfEachNode() const;

    const int getDegreeCount() const;

    const int getEdgeCount() const;

    const int getNodeCount() const;

    const std::vector<std::vector<int>>& getAdjacencyList() const;

    const std::map<int, int>& getMapping() const;

private:

    void readAdjacencyList( const std::string& fileName );

    void sortNeighboursInAscendingOrderForEachNode( std::vector<std::vector<int>>& net );

    void setStubs( std::vector<int>& stubs );

    void shuffleStubs( std::vector<int>& stubs );

    std::vector<double> normaliseDegreeDistribution( std::vector<int>& distribution );

    double calculateDegreeDistributionMean( const std::vector<int>& distribution );

    double getAverageExcessDegree() const;

    double calculateExcessDegreeVariance( double mean ) const;

    double getClosenessCentrality( std::vector<int> distance ) const;

    std::vector<int> getNeighboursShuffledThenSortedByDegreeByDesc( int current ) const;

    int getTargetOrRandomHighestDegreeWhiteOrGreyNeighbour( int current,
                                                            int target,
                                                            std::vector<Colour>& colours ) const;

    void setNodeColour( std::vector<Colour>& colours, int current ) const;

    bool isCurrentNodeNeighbourOfTarget( int current, int target ) const;

    int getRandomHighestDegreeNeighbour( int current );

    std::vector<int> getNodes() const;

    std::vector<int> getNeighboursByDegree( int current ) const;

    std::vector<int> getLayer( const std::vector<int>& nodes,
                               const std::vector<int>& degreeOfNeighbours,
                               int core ) const;

    int getMinDegree( std::vector<int> degreeOfNodes ) const;

    std::vector<std::vector<int>> adjacencyList;
    int degreeCount;
    int edgeCount;
    int nodeCount;

    std::map<int, int> mapping;
};

#endif //ANALYSEADJACENCYLIST_H
