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
package com.graphhopper;

import com.graphhopper.config.CHProfile;
import com.graphhopper.config.Profile;
import com.graphhopper.routing.MultiplePointsNotFoundException;
import com.graphhopper.routing.matrix.DistanceMatrix;
import com.graphhopper.routing.matrix.GHMatrixRequest;
import com.graphhopper.util.Helper;
import com.graphhopper.util.shapes.GHPoint;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvFileSource;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static com.graphhopper.routing.matrix.DistanceMatrix.DISTANCE_SNAP_ERROR_VALUE;
import static com.graphhopper.routing.matrix.DistanceMatrix.TIME_SNAP_ERROR_VALUE;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class GraphHopperMatrixTest {
    private static final double DISTANCE_METERS_DELTA = 1d;
    private static final double TIME_MILLIS_DELTA = 10d;

    public static final String DIR = "../core/files";

    // map locations
    private static final String ANDORRA = DIR + "/andorra.osm.gz";
    private static final String MONACO = DIR + "/monaco.osm.gz";

    private static final String NOTHINGHAM = DIR + "/nottinghamshire.osm.pbf";

    // when creating GH instances make sure to use this as the GH location such that it will be cleaned between tests
    private static final String GH_LOCATION = "target/graphhopper-test-gh";

    @BeforeEach
    @AfterEach
    public void setup() {
        Helper.removeDir(new File(GH_LOCATION));
    }

    static class MatrixText {
        private List<GHPoint> origins = new ArrayList<>();
        private List<GHPoint> destinations = new ArrayList<>();

        public List<GHPoint> getOrigins() {
            return origins;
        }

        public void setOrigins(List<GHPoint> origins) {
            this.origins = origins;
        }

        public List<GHPoint> getDestinations() {
            return destinations;
        }

        public void setDestinations(List<GHPoint> destinations) {
            this.destinations = destinations;
        }
    }

    private MatrixText parseRawMatrix(String rawMatrix) {
        String[] rawPoints = rawMatrix.split("&");
        Stream<GHPoint> sourceGHPoints = Arrays.stream(rawPoints).filter(p -> p.startsWith("s")).map(p -> {
            String[] latLong = p.replace("s=", "").split(",");
            return new GHPoint(Double.parseDouble(latLong[0]), Double.parseDouble(latLong[1]));
        });
        Stream<GHPoint> destGHPoints = Arrays.stream(rawPoints).filter(p -> p.startsWith("d")).map(p -> {
            String[] latLong = p.replace("d=", "").split(",");
            return new GHPoint(Double.parseDouble(latLong[0]), Double.parseDouble(latLong[1]));
        });

        MatrixText test = new MatrixText();
        test.setDestinations(destGHPoints.collect(Collectors.toList()));
        test.setOrigins(sourceGHPoints.collect(Collectors.toList()));
        return test;
    }

    @ParameterizedTest(name = "Andorra Matrix Test {index} => matrix {1}")
    @CsvFileSource(resources = "/com/graphhopper/routing/matrix/andorra_matrix.csv", delimiter = ';')
    public void testAndorraDistanceMatrixResultsConsistencyAgainstRoutePoint2Point(String country, String matrix) {
        testDistanceMatrixResultsConsistency(matrix, ANDORRA);
    }

    @ParameterizedTest(name = "Montecarlo Matrix Test {index} => matrix {1}")
    @CsvFileSource(resources = "/com/graphhopper/routing/matrix/monaco_matrix.csv", delimiter = ';')
    public void testMonacoDistanceMatrixResultsConsistencyAgainstRoutePoint2Point(String country, String matrix) {
        testDistanceMatrixResultsConsistency(matrix, MONACO);
    }

    private void testDistanceMatrixResultsConsistency(String matrix, String osmLocation) {
        final MatrixText matrixText = parseRawMatrix(matrix);

        Profile carProfile = new Profile("car");
        carProfile.setTurnCosts(false);
        CHProfile chCarProfile = new CHProfile("car");

        List<Profile> profiles = new ArrayList<>();
        profiles.add(carProfile);

        List<CHProfile> chProfiles = new ArrayList<>();
        chProfiles.add(chCarProfile);

        GraphHopperConfig config = new GraphHopperConfig();
        config.setProfiles(profiles);
        config.setCHProfiles(chProfiles);


        GraphHopper hopper = new GraphHopper()
                .setOSMFile(osmLocation)
                .init(config)
                .setGraphHopperLocation(GH_LOCATION)
                .importOrLoad();

        List<GHPoint> origins = matrixText.getOrigins();
        List<GHPoint> destinations = matrixText.getDestinations();

        GHMatrixRequest request = new GHMatrixRequest();
        request.setProfile("car");
        request.setOrigins(origins);
        request.setDestinations(destinations);

        DistanceMatrix matrixResult = hopper.matrix(request).getMatrix();

        for (int sourceIdx = 0; sourceIdx < origins.size(); sourceIdx++) {
            for (int destIdx = 0; destIdx < destinations.size(); destIdx++) {
                double dist = matrixResult.getDistance(sourceIdx, destIdx);
                double time = matrixResult.getTime(sourceIdx, destIdx);

                List<GHPoint> points = new ArrayList<>();
                points.add(origins.get(sourceIdx));
                points.add(destinations.get(destIdx));

                GHRequest req = new GHRequest().setProfile("car").setPoints(points);
                GHResponse response = hopper.route(req);

                assertEquals(response.getBest().getDistance(), dist, DISTANCE_METERS_DELTA, String.format("Unexpected distance for points=%s", points));
                assertEquals(response.getBest().getTime(), time, TIME_MILLIS_DELTA, String.format("Unexpected time for points=%s", points));
            }
        }
    }

    @Test
    void testFailFast() {
        List<GHPoint> originPoints = new ArrayList<>();
        originPoints.add(new GHPoint(42.50488142419136, 1.5239627305424333));
        originPoints.add(new GHPoint(42.500908692810256, 1.0666521515059635));
        List<GHPoint> destinationPoints = new ArrayList<>();
        destinationPoints.add(new GHPoint(42.50272031776507, 1.5148522642488023));
        destinationPoints.add(new GHPoint(42.50272031776507, 1.5148522642488023));

        final MatrixText matrixText = new MatrixText();
        matrixText.setDestinations(originPoints);
        matrixText.setOrigins(destinationPoints);

        Profile carProfile = new Profile("car");
        carProfile.setTurnCosts(false);
        CHProfile chCarProfile = new CHProfile("car");

        List<Profile> profiles = new ArrayList<>();
        profiles.add(carProfile);

        List<CHProfile> chProfiles = new ArrayList<>();
        chProfiles.add(chCarProfile);

        GraphHopperConfig config = new GraphHopperConfig();
        config.setProfiles(profiles);
        config.setCHProfiles(chProfiles);


        GraphHopper hopper = new GraphHopper()
                .setOSMFile(ANDORRA)
                .init(config)
                .setGraphHopperLocation(GH_LOCATION)
                .importOrLoad();

        List<GHPoint> origins = matrixText.getOrigins();
        List<GHPoint> destinations = matrixText.getDestinations();

        GHMatrixRequest request = new GHMatrixRequest();
        request.setProfile("car");
        request.setOrigins(origins);
        request.setDestinations(destinations);

        MultiplePointsNotFoundException thrown = Assertions.assertThrows(MultiplePointsNotFoundException.class, () -> {
            hopper.matrix(request).getMatrix();
        });

        Assertions.assertTrue(thrown.getPointsNotFound().size() > 0);
    }

    @Test
    void testFailFastDisabled() {
        List<GHPoint> originPoints = new ArrayList<>();
        originPoints.add(new GHPoint(42.50488142419136, 1.5239627305424333));
        originPoints.add(new GHPoint(42.500908692810256, 1.0666521515059635)); // Point outside Andorra (Not found)
        originPoints.add(new GHPoint(42.50488142419136, 1.5239627305424333));
        originPoints.add(new GHPoint(42.500908692810256, 1.0666521515059635)); // Point outside Andorra (Not found)
        List<GHPoint> destinationPoints = new ArrayList<>();
        destinationPoints.add(new GHPoint(42.50272031776507, 1.5148522642488023));
        destinationPoints.add(new GHPoint(42.50272031776507, 1.5148522642488023));

        final MatrixText matrixText = new MatrixText();
        matrixText.setDestinations(originPoints);
        matrixText.setOrigins(destinationPoints);

        Profile carProfile = new Profile("car");
        carProfile.setTurnCosts(false);
        CHProfile chCarProfile = new CHProfile("car");

        List<Profile> profiles = new ArrayList<>();
        profiles.add(carProfile);

        List<CHProfile> chProfiles = new ArrayList<>();
        chProfiles.add(chCarProfile);

        GraphHopperConfig config = new GraphHopperConfig();
        config.setProfiles(profiles);
        config.setCHProfiles(chProfiles);


        GraphHopper hopper = new GraphHopper()
                .setOSMFile(ANDORRA)
                .init(config)
                .setGraphHopperLocation(GH_LOCATION)
                .importOrLoad();

        List<GHPoint> origins = matrixText.getOrigins();
        List<GHPoint> destinations = matrixText.getDestinations();

        GHMatrixRequest request = new GHMatrixRequest();
        request.setProfile("car");
        request.setFailFast(false);
        request.setOrigins(origins);
        request.setDestinations(destinations);

        DistanceMatrix result = hopper.matrix(request).getMatrix();
        System.out.println(result);
        // Distance assertions
        assertEquals(DISTANCE_SNAP_ERROR_VALUE, result.getDistance(0, 1));
        assertEquals(DISTANCE_SNAP_ERROR_VALUE, result.getDistance(1, 1));
        assertTrue(result.getDistance(0, 0) > 0);
        assertTrue(result.getDistance(1, 0) > 0);
        // Time assertions
        assertEquals(TIME_SNAP_ERROR_VALUE, result.getTime(0, 1));
        assertEquals(TIME_SNAP_ERROR_VALUE, result.getTime(1, 1));
        assertTrue(result.getTime(0, 0) > 0);
        assertTrue(result.getTime(1, 0) > 0);
    }

    @Test
    void testSpeedZeroEdge(){

        List<GHPoint> originPoints = new ArrayList<>();
        originPoints.add(new GHPoint(52.95675659,-1.14968896));

        List<GHPoint> destinationPoints = new ArrayList<>();
        destinationPoints.add(new GHPoint(52.98455,-1.19811));

        final MatrixText matrixText = new MatrixText();
        matrixText.setDestinations(originPoints);
        matrixText.setOrigins(destinationPoints);

        Profile carProfile = new Profile("car");
        carProfile.setTurnCosts(false);
        CHProfile chCarProfile = new CHProfile("car");

        List<Profile> profiles = new ArrayList<>();
        profiles.add(carProfile);

        List<CHProfile> chProfiles = new ArrayList<>();
        chProfiles.add(chCarProfile);

        GraphHopperConfig config = new GraphHopperConfig();
        config.setProfiles(profiles);
        config.setCHProfiles(chProfiles);


        GraphHopper hopper = new GraphHopper()
                .setOSMFile(NOTHINGHAM)
                .init(config)
                .setGraphHopperLocation(GH_LOCATION)
                .importOrLoad();

        List<GHPoint> origins = matrixText.getOrigins();
        List<GHPoint> destinations = matrixText.getDestinations();

        GHMatrixRequest request = new GHMatrixRequest();
        request.setProfile("car");
        request.setFailFast(false);
        request.setOrigins(origins);
        request.setDestinations(destinations);

        DistanceMatrix result = hopper.matrix(request).getMatrix();
        System.out.println(result);
        // Distance assertions
        assertTrue(Math.round(result.getDistance(0, 0)) == 5260);
        // Time assertions
        assertTrue(result.getTime(0, 0) == 414191);

    }
}
