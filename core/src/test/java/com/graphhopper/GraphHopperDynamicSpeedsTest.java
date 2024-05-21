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
import com.graphhopper.config.LMProfile;
import com.graphhopper.config.Profile;
import com.graphhopper.reader.ReaderWay;
import com.graphhopper.reader.dem.SRTMProvider;
import com.graphhopper.reader.dem.SkadiProvider;
import com.graphhopper.routing.ev.EdgeIntAccess;
import com.graphhopper.routing.ev.EncodedValueLookup;
import com.graphhopper.routing.ev.RoadEnvironment;
import com.graphhopper.routing.ev.Subnetwork;
import com.graphhopper.routing.util.AllEdgesIterator;
import com.graphhopper.routing.util.DefaultSnapFilter;
import com.graphhopper.routing.util.EdgeFilter;
import com.graphhopper.routing.util.countryrules.CountryRuleFactory;
import com.graphhopper.routing.util.parsers.DefaultTagParserFactory;
import com.graphhopper.routing.util.parsers.OSMRoadEnvironmentParser;
import com.graphhopper.routing.util.parsers.TagParser;
import com.graphhopper.routing.weighting.Weighting;
import com.graphhopper.search.KVStorage;
import com.graphhopper.speeds.WaySpeedsProvider;
import com.graphhopper.storage.IntsRef;
import com.graphhopper.storage.index.LocationIndexTree;
import com.graphhopper.storage.index.Snap;
import com.graphhopper.util.*;
import com.graphhopper.util.Parameters.CH;
import com.graphhopper.util.Parameters.Landmark;
import com.graphhopper.util.Parameters.Routing;
import com.graphhopper.util.details.PathDetail;
import com.graphhopper.util.exceptions.ConnectionNotFoundException;
import com.graphhopper.util.exceptions.MaximumNodesExceededException;
import com.graphhopper.util.exceptions.PointDistanceExceededException;
import com.graphhopper.util.shapes.BBox;
import com.graphhopper.util.shapes.GHPoint;
import com.graphhopper.util.shapes.GHPoint3D;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.ValueSource;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;

import java.io.File;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

import static com.graphhopper.json.Statement.If;
import static com.graphhopper.json.Statement.Op.MULTIPLY;
import static com.graphhopper.util.GHUtility.createCircle;
import static com.graphhopper.util.GHUtility.createRectangle;
import static com.graphhopper.util.Parameters.Algorithms.*;
import static com.graphhopper.util.Parameters.Curbsides.*;
import static com.graphhopper.util.Parameters.Routing.TIMEOUT_MS;
import static com.graphhopper.util.Parameters.Routing.U_TURN_COSTS;
import static java.util.Arrays.asList;
import static org.junit.jupiter.api.Assertions.*;

/**
 * @author Peter Karich
 */
public class GraphHopperDynamicSpeedsTest {

    public static final String DIR = "../core/files";

    // map locations
    private static final String MONACO = DIR + "/monaco.osm.gz";

    // when creating GH instances make sure to use this as the GH location such that it will be cleaned between tests
    private static final String GH_LOCATION = "target/graphhopper-test-gh";
    private static final String GH_LOCATION_CUSTOM_SPEEDS = "target/graphhopper-test-gh-custom_speeds";

    @BeforeEach
    @AfterEach
    public void setup() {

        Helper.removeDir(new File(GH_LOCATION));
        Helper.removeDir(new File(GH_LOCATION_CUSTOM_SPEEDS));
    }

    class CustomWaySpeedProvider implements WaySpeedsProvider {

        @Override
        public double speedForWay(long osmWayId) {
            return 10.1;
        }
    }

    @Test
    public void testMultipleVehiclesWithCH() {
        final String bikeProfile = "bike_profile";
        final String carProfile = "car_profile";
        List<Profile> profiles = asList(
                new Profile(bikeProfile).setVehicle("bike"),
                new Profile(carProfile).setVehicle("car")
        );

        GraphHopper hopper = new GraphHopper().
                setGraphHopperLocation(GH_LOCATION).
                setOSMFile(MONACO).
                setProfiles(profiles).
                setStoreOnFlush(true);
        hopper.getCHPreparationHandler().setCHProfiles(
                new CHProfile(bikeProfile),
                new CHProfile(carProfile)
        );

        GraphHopper hopperCustomSpeeds = new GraphHopper().
                setGraphHopperLocation(GH_LOCATION_CUSTOM_SPEEDS).
                setOSMFile(MONACO).
                setProfiles(profiles).
                setStoreOnFlush(true)
                .setEncodedValuesString("roundabout, road_class, road_class_link, road_environment, max_speed, road_access, ferry_speed, bike_network, get_off_bike, smoothness, osm_way_id");
        hopperCustomSpeeds.getCHPreparationHandler().setCHProfiles(
                new CHProfile(bikeProfile),
                new CHProfile(carProfile)
        );

        hopperCustomSpeeds.setDynamicSpeeds(new CustomWaySpeedProvider());


        hopper.importOrLoad();
        hopperCustomSpeeds.importOrLoad();

        GHResponse rsp = hopper.route(new GHRequest(43.73005, 7.415707, 43.741522, 7.42826)
                .setProfile(carProfile));
        ResponsePath res = rsp.getBest();
        assertFalse(rsp.hasErrors(), rsp.getErrors().toString());
        assertEquals(205, res.getTime() / 1000f, 1);
        assertEquals(2837, res.getDistance(), 1);

        GHResponse rspSpeeds = hopperCustomSpeeds.route(new GHRequest(43.73005, 7.415707, 43.741522, 7.42826)
                .setProfile(carProfile));
        ResponsePath resSpeeds = rspSpeeds.getBest();
        assertFalse(rspSpeeds.hasErrors(), rspSpeeds.getErrors().toString());
        assertEquals(893, resSpeeds.getTime() / 1000f, 1);
        assertEquals(2481, resSpeeds.getDistance(), 1);


        GHResponse rspBike = hopper.route(new GHRequest(43.73005, 7.415707, 43.741522, 7.42826)
                .setProfile(bikeProfile));
        ResponsePath resBike = rspBike.getBest();
        assertFalse(rspBike.hasErrors(), rspBike.getErrors().toString());
        assertEquals(536, resBike.getTime() / 1000f, 1);
        assertEquals(2521, resBike.getDistance(), 1);

        GHResponse rspBikeSpeeds = hopperCustomSpeeds.route(new GHRequest(43.73005, 7.415707, 43.741522, 7.42826)
                .setProfile(bikeProfile));
        ResponsePath resBikeSpeeds = rspBikeSpeeds.getBest();
        assertFalse(rspBikeSpeeds.hasErrors(), rspBikeSpeeds.getErrors().toString());
        assertEquals(834, resBikeSpeeds.getTime() / 1000f, 1);
        assertEquals(2318, resBikeSpeeds.getDistance(), 1);
        
    }

}

