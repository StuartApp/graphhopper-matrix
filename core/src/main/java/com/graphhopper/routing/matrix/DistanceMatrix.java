package com.graphhopper.routing.matrix;


import com.carrotsearch.hppc.IntArrayList;

import java.util.*;
import java.util.function.Consumer;

/**
 * Holds the resulting distance matrix for a given set of
 * origin/destination nodes.
 * <p>
 * The matrix consists of origin points which form the rows,
 * and destination points which form the columns.
 *
 * @author Pascal Büttiker
 */
public final class DistanceMatrix {

    public static final double DISTANCE_SNAP_ERROR_VALUE = Double.MAX_VALUE;
    public static final long TIME_SNAP_ERROR_VALUE = Long.MAX_VALUE;

    private final int numberOfOrigins;
    private final int numberOfDestinations;

    private final double[][] distances;

    private final long[][] times;

    private IntArrayList snapOriginsErrors;
    private IntArrayList snapDestinationErrors;

    public double[][] getDistances() {
        return distances;
    }

    public long[][] getTimes() {
        return times;
    }

    /**
     * Creates a new GHMatrixResponse with the given dimensions
     *
     * @param numberOfOrigins      The number of origin points (rows)
     * @param numberOfDestinations The number of destination points (columns)
     */

    public DistanceMatrix(int numberOfOrigins, int numberOfDestinations, IntArrayList originsErrors, IntArrayList destinationsErros) {
        this.numberOfOrigins = numberOfOrigins;
        this.numberOfDestinations = numberOfDestinations;

        distances = new double[numberOfOrigins][numberOfDestinations];
        times = new long[numberOfOrigins][numberOfDestinations];

        this.snapOriginsErrors = originsErrors;
        this.snapDestinationErrors = destinationsErros;

        for(int i = 0; i < numberOfOrigins; i++){
            setSnapErrorsValues(i);
        }

    }

    private void setSnapErrorValue(int origin, int destination){
        distances[origin][destination] = DISTANCE_SNAP_ERROR_VALUE;
        times[origin][destination] = TIME_SNAP_ERROR_VALUE;
    }

    private void setSnapErrorsValues(int originIdx){

        boolean originError = snapOriginsErrors.contains(originIdx);
        for(int i = 0; i < numberOfDestinations; i++){
            if(originError || snapDestinationErrors.contains(i)){
                setSnapErrorValue(originIdx,i);
            }
        }
    }

    /**
     * Gets the number of origins
     * (rows)
     */
    public int getNumberOfOrigins() {
        return numberOfOrigins;
    }

    /**
     * Gets the number of destinations
     * (columns)
     */
    public int getNumberOfDestinations() {
        return numberOfDestinations;
    }

    /**
     * Set the distance/time info of a single (origin --> destination) cell
     *
     * @param originIndex The index of the origin
     * @param destIndex   The index of the destination
     * @param distance    Distance value
     * @param time        Time value
     */
    public void setCell(int originIndex, int destIndex, double distance, long time) {
        distances[originIndex][destIndex] = distance;
        times[originIndex][destIndex] = time;
    }


    public double getDistance(int originIndex, int destIndex) {
        return distances[originIndex][destIndex];
    }

    public long getTime(int originIndex, int destIndex) {
        return times[originIndex][destIndex];
    }

    /**
     * Returns the matrix as formatted string
     */
    @Override
    public String toString() {
        String allMatrices = "";
        if (distances != null) {
            allMatrices += Arrays.deepToString(distances) + "\n";
        }
        if (times != null) {
            allMatrices += Arrays.deepToString(times) + "\n";
        }
        return allMatrices;
    }
}
