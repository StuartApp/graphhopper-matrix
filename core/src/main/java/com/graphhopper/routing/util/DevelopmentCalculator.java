/*
 *  Licensed to GraphHopper GmbH under one or more contributor
 *  license agreements. See the NOTICE file distributed with this work for
 *  additional information regarding copyright ownership.
 *
 *  GraphHopper GmbH licenses this file to you under the Apache License,
 *  Version 2.0 (the "License"); you may not use this file except in
 *  compliance with the License. You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

package com.graphhopper.routing.util;

import com.graphhopper.routing.ev.BooleanEncodedValue;
import com.graphhopper.routing.ev.Development;
import com.graphhopper.routing.ev.EnumEncodedValue;
import com.graphhopper.routing.ev.RoadClass;
import com.graphhopper.storage.Graph;
import com.graphhopper.util.EdgeIteratorState;
import com.graphhopper.util.StopWatch;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.function.ToDoubleFunction;

public class DevelopmentCalculator {
    private static final Logger logger = LoggerFactory.getLogger(RoadDensityCalculator.class);

    /**
     * Calculates the development (rural/residential/city) for all edges of the graph.
     * First a weighted road density is calculated for every edge to determine whether it belongs to a residential area.
     * In a second step very dense residential areas are classified as 'city'.
     *
     * @param residentialAreaRadius      radius used for residential area calculation in meters
     * @param residentialAreaSensitivity Use higher values if there are too many residential areas that are not recognized. Use
     *                                   smaller values if there are too many non-residential areas that are classified as residential.
     * @param cityAreaRadius             in meters, see residentialAreaRadius
     * @param cityAreaSensitivity        similar to residentialAreaSensitivity, but for the city classification
     * @param threads                    number of threads used to calculate the road densities
     */
    public static void calcDevelopment(Graph graph, EnumEncodedValue<Development> developmentEnc,
                                       EnumEncodedValue<RoadClass> roadClassEnc, BooleanEncodedValue roadClassLinkEnc,
                                       double residentialAreaRadius, double residentialAreaSensitivity,
                                       double cityAreaRadius, double cityAreaSensitivity,
                                       int threads) {
        logger.info("Calculating residential areas ..., radius={}, sensitivity={}, threads={}", residentialAreaRadius, residentialAreaSensitivity, threads);
        StopWatch sw = StopWatch.started();
        calcResidential(graph, developmentEnc, roadClassEnc, roadClassLinkEnc, residentialAreaRadius, residentialAreaSensitivity, threads);
        logger.info("Finished calculating residential areas, took: " + sw.stop().getSeconds() + "s");
        if (cityAreaRadius > 1) {
            logger.info("Calculating city areas ..., radius={}, sensitivity={}, threads={}", cityAreaRadius, cityAreaSensitivity, threads);
            sw = StopWatch.started();
            calcCity(graph, developmentEnc, cityAreaRadius, cityAreaSensitivity, threads);
            logger.info("Finished calculating city areas, took: " + sw.stop().getSeconds() + "s");
        }
    }

    private static void calcResidential(Graph graph, EnumEncodedValue<Development> developmentEnc,
                                        EnumEncodedValue<RoadClass> roadClassEnc, BooleanEncodedValue roadClassLinkEnc,
                                        double radius, double sensitivity, int threads) {
        final ToDoubleFunction<EdgeIteratorState> calcRoadFactor = edge -> {
            if (edge.get(roadClassLinkEnc))
                // highway exits sometimes have a 'high' density. By excluding them here we try to prevent them from
                // being classified as residential areas when they are in the countryside.
                return 0;
            RoadClass roadClass = edge.get(roadClassEnc);
            switch (roadClass) {
                // we're interested in the road density of urban roads, so residential areas are particularly interesting
                case RESIDENTIAL:
                case LIVING_STREET:
                case FOOTWAY:
                case CYCLEWAY:
                case STEPS:
                    return 2;
                case MOTORWAY:
                case TRUNK:
                case PRIMARY:
                case SECONDARY:
                case TERTIARY:
                case SERVICE:
                case OTHER:
                    return 1;
                default:
                    return 0;
            }
        };
        // temporarily write results to an external array for thread-safety
        boolean[] isResidential = new boolean[graph.getEdges()];
        RoadDensityCalculator.calcRoadDensities(graph, (calculator, edge) -> {
            RoadClass roadClass = edge.get(roadClassEnc);
            if (roadClass == RoadClass.RESIDENTIAL || roadClass == RoadClass.LIVING_STREET) {
                isResidential[edge.getEdge()] = true;
                return;
            }
            double roadDensity = calculator.calcRoadDensity(edge, radius, calcRoadFactor);
            isResidential[edge.getEdge()] = roadDensity * sensitivity >= 1.0;
        }, threads);
        for (int edge = 0; edge < isResidential.length; edge++)
            graph.getEdgeIteratorState(edge, Integer.MIN_VALUE).set(developmentEnc, isResidential[edge] ? Development.RESIDENTIAL : Development.RURAL);
    }

    private static void calcCity(Graph graph, EnumEncodedValue<Development> developmentEnc,
                                 double radius, double sensitivity, int threads) {
        // do not modify the development values as long as we are still reading them -> store city flags in this array first
        boolean[] isCity = new boolean[graph.getEdges()];
        final ToDoubleFunction<EdgeIteratorState> calcRoadFactor = edge -> edge.get(developmentEnc) == Development.RESIDENTIAL ? 1 : 0;
        RoadDensityCalculator.calcRoadDensities(graph, (calculator, edge) -> {
            Development development = edge.get(developmentEnc);
            if (development == Development.RURAL)
                return;
            double roadDensity = calculator.calcRoadDensity(edge, radius, calcRoadFactor);
            if (roadDensity * sensitivity >= 1.0)
                isCity[edge.getEdge()] = true;
        }, threads);
        for (int edge = 0; edge < isCity.length; edge++)
            if (isCity[edge])
                graph.getEdgeIteratorState(edge, Integer.MIN_VALUE).set(developmentEnc, Development.CITY);
    }
}
