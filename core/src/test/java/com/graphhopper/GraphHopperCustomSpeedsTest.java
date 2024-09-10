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
import com.graphhopper.routing.ev.RoadClass;
import com.graphhopper.routing.ev.VehicleSpeed;
import com.graphhopper.speeds.SpeedKmByHour;
import com.graphhopper.speeds.WaySpeedsProvider;
import com.graphhopper.util.*;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.util.*;

import static java.util.Arrays.asList;
import static org.junit.jupiter.api.Assertions.*;

/**
 * @author Peter Karich
 */
public class GraphHopperCustomSpeedsTest {

    public static final String DIR = "../core/files";

    // map locations
    private static final String MONACO = DIR + "/monaco.osm.gz";

    // when creating GH instances make sure to use this as the GH location such that it will be cleaned between tests
    private static final String GH_LOCATION = "target/graphhopper-test-gh";
    private static final String GH_LOCATION_CUSTOM_SPEEDS = "target/graphhopper-test-gh-custom_speeds";

    private static final String GH_LOCATION_CUSTOM_SPEEDS_RELOAD = "target/graphhopper-test-gh-custom_speeds-reload";

    @BeforeEach
    @AfterEach
    public void setup() {

        Helper.removeDir(new File(GH_LOCATION));
        Helper.removeDir(new File(GH_LOCATION_CUSTOM_SPEEDS));
    }

    @Test
    public void testDynamicSpeedProvider() {

        class CustomWaySpeedProvider implements WaySpeedsProvider {

            @Override
            public Optional<SpeedKmByHour> speedForWay(long osmWayId) {
                return Optional.of(new SpeedKmByHour(10.10));
            }

            @Override
            public Optional<SpeedKmByHour> speedForRoadClass(RoadClass roadClass) {
                return Optional.of(new SpeedKmByHour(10.1));
            }
        }

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

        GraphHopper hopperCustomSpeeds = new GraphHopperCustomSpeeds(new CustomWaySpeedProvider()).
                setGraphHopperLocation(GH_LOCATION_CUSTOM_SPEEDS).
                setOSMFile(MONACO).
                setProfiles(profiles).
                setStoreOnFlush(true)
                .setEncodedValuesString("roundabout, road_class, road_class_link, road_environment, max_speed, road_access, ferry_speed, bike_network, get_off_bike, smoothness, osm_way_id");
        hopperCustomSpeeds.getCHPreparationHandler().setCHProfiles(
                new CHProfile(bikeProfile),
                new CHProfile(carProfile)
        );

        hopper.importOrLoad();
        hopperCustomSpeeds.importOrLoad();

        GHRequest request = new GHRequest(43.73005, 7.415707, 43.741522, 7.42826);

        GHResponse rsp = hopper.route(request.setProfile(carProfile));
        ResponsePath res = rsp.getBest();
        assertFalse(rsp.hasErrors(), rsp.getErrors().toString());
        assertEquals(205, res.getTime() / 1000f, 1);
        assertEquals(2837, res.getDistance(), 1);

        GHResponse rspSpeeds = hopperCustomSpeeds.route(request.setProfile(carProfile));
        ResponsePath resSpeeds = rspSpeeds.getBest();
        assertFalse(rspSpeeds.hasErrors(), rspSpeeds.getErrors().toString());
        assertEquals(893, resSpeeds.getTime() / 1000f, 1);
        assertEquals(2481, resSpeeds.getDistance(), 1);


        GHResponse rspBike = hopper.route(request.setProfile(bikeProfile));
        ResponsePath resBike = rspBike.getBest();
        assertFalse(rspBike.hasErrors(), rspBike.getErrors().toString());
        assertEquals(536, resBike.getTime() / 1000f, 1);
        assertEquals(2521, resBike.getDistance(), 1);

        GHResponse rspBikeSpeeds = hopperCustomSpeeds.route(request.setProfile(bikeProfile));
        ResponsePath resBikeSpeeds = rspBikeSpeeds.getBest();
        assertFalse(rspBikeSpeeds.hasErrors(), rspBikeSpeeds.getErrors().toString());
        assertEquals(834, resBikeSpeeds.getTime() / 1000f, 1);
        assertEquals(2318, resBikeSpeeds.getDistance(), 1);

    }

    @Test
    public void testDynamicSpeedProviderDiscardCustomSpeedIfHigherThanMaxSpeed() {

        class SpeedsTooHighSpeedProvider implements WaySpeedsProvider {

            @Override
            public Optional<SpeedKmByHour> speedForWay(long osmWayId) {
                return Optional.of(new SpeedKmByHour(100));
            }

            @Override
            public Optional<SpeedKmByHour> speedForRoadClass(RoadClass roadClass) {
                return Optional.of(new SpeedKmByHour(80.1));
            }
        }

        final String carProfile = "car_profile";

        List<Profile> profiles = asList(new Profile(carProfile).setVehicle("car"));

        GraphHopper hopper = new GraphHopper().
                setGraphHopperLocation(GH_LOCATION).
                setOSMFile(MONACO).
                setProfiles(profiles).
                setStoreOnFlush(true);
        hopper.getCHPreparationHandler().setCHProfiles(
                new CHProfile(carProfile)
        );

        GraphHopper hopperCustomSpeeds = new GraphHopperCustomSpeeds(new SpeedsTooHighSpeedProvider()).
                setGraphHopperLocation(GH_LOCATION_CUSTOM_SPEEDS).
                setOSMFile(MONACO).
                setProfiles(profiles).
                setStoreOnFlush(true)
                .setEncodedValuesString("roundabout, road_class, road_class_link, road_environment, max_speed, road_access, ferry_speed, bike_network, get_off_bike, smoothness, osm_way_id");
        hopperCustomSpeeds.getCHPreparationHandler().setCHProfiles(
                new CHProfile(carProfile)
        );

        hopper.importOrLoad();
        hopperCustomSpeeds.importOrLoad();

        // requests uses osm way 166009792 (secondary road, that has max speed = 50km)
        GHRequest request = new GHRequest(43.73887, 7.42074, 43.73992, 7.42386);

        GHResponse rsp = hopper.route(request.setProfile(carProfile));
        ResponsePath res = rsp.getBest();
        assertFalse(rsp.hasErrors(), rsp.getErrors().toString());
        assertEquals(21, res.getTime() / 1000f, 1); // 47kmh = 13.2m/s
        assertEquals(277, res.getDistance(), 1);

        GHResponse rspSpeeds = hopperCustomSpeeds.route(request
                .setProfile(carProfile));
        ResponsePath resSpeeds = rspSpeeds.getBest();
        assertFalse(rspSpeeds.hasErrors(), rspSpeeds.getErrors().toString());
        assertEquals(21, resSpeeds.getTime() / 1000f, 1);  // 47kmh = 13.2m/s
        assertEquals(277, resSpeeds.getDistance(), 1);

        // the result with and without custom speeds should be the same when the custom speed is > than the max_speed
        assertEquals(res.getTime() / 1000f, resSpeeds.getTime() / 1000f, 1);


    }

    @Test
    public void testDynamicSpeedProviderDiscardCustomSpeedIfHigherThanMaxEncoderSpeed() {

        // For bike, the max encoder speed is 30kmh
        // For car, the max encoder speed is 254 kmh
        // We are trying to set a custom speed for bike > 30kmh, that is the limit for the bike
        // encoder, so we expect it to be discarded, and we expect that the default speed
        // for those roads will be used

        class SpeedsTooHighSpeedProvider implements WaySpeedsProvider {

            @Override
            public Optional<SpeedKmByHour> speedForWay(long osmWayId) {
                return Optional.of(new SpeedKmByHour(100));
            }

            @Override
            public Optional<SpeedKmByHour> speedForRoadClass(RoadClass roadClass) {
                return Optional.of(new SpeedKmByHour(80.1));
            }
        }

        final String bikeProfile = "bike_profile";

        List<Profile> profiles = asList(
                new Profile(bikeProfile).setVehicle("bike")
        );

        GraphHopper hopper = new GraphHopper().
                setGraphHopperLocation(GH_LOCATION).
                setOSMFile(MONACO).
                setProfiles(profiles).
                setStoreOnFlush(true);
        hopper.getCHPreparationHandler().setCHProfiles(
                new CHProfile(bikeProfile)
        );

        GraphHopper hopperCustomSpeeds = new GraphHopperCustomSpeeds(new SpeedsTooHighSpeedProvider()).
                setGraphHopperLocation(GH_LOCATION_CUSTOM_SPEEDS).
                setOSMFile(MONACO).
                setProfiles(profiles).
                setStoreOnFlush(true)
                .setEncodedValuesString("roundabout, road_class, road_class_link, road_environment, max_speed, road_access, ferry_speed, bike_network, get_off_bike, smoothness, osm_way_id");
        hopperCustomSpeeds.getCHPreparationHandler().setCHProfiles(
                new CHProfile(bikeProfile)
        );

        hopper.importOrLoad();
        hopperCustomSpeeds.importOrLoad();

        GHRequest request = new GHRequest(43.73887, 7.42074, 43.73992, 7.42386);

        GHResponse rspBike = hopper.route(request
                .setProfile(bikeProfile));
        ResponsePath resBike = rspBike.getBest();
        assertFalse(rspBike.hasErrors(), rspBike.getErrors().toString());
        assertEquals(55, resBike.getTime() / 1000f, 1); // 18kmh = 5 m/s
        assertEquals(277, resBike.getDistance(), 1);

        GHResponse rspBikeSpeeds = hopperCustomSpeeds.route(request
                .setProfile(bikeProfile));
        ResponsePath resBikeSpeeds = rspBikeSpeeds.getBest();
        assertFalse(rspBikeSpeeds.hasErrors(), rspBikeSpeeds.getErrors().toString());
        assertEquals(55, resBikeSpeeds.getTime() / 1000f, 1); // 17kmh = 4.7 m/s
        assertEquals(277, resBikeSpeeds.getDistance(), 1);

        // the result with and without custom speeds should be the same when the custom speed is > than the encoder limit
        assertEquals(resBike.getTime() / 1000f, resBikeSpeeds.getTime() / 1000f, 1);

    }

    @Test
    public void testEncodingManagerStorePropertiesAfterSetCustomSpeeds() {

        // For bike, the max encoder speed is 30kmh
        class SpeedsTooHighSpeedProvider implements WaySpeedsProvider {

            @Override
            public Optional<SpeedKmByHour> speedForWay(long osmWayId) {
                return Optional.of(new SpeedKmByHour(28));
            }

            @Override
            public Optional<SpeedKmByHour> speedForRoadClass(RoadClass roadClass) {
                return Optional.of(new SpeedKmByHour(28));
            }
        }

        final String bikeProfile = "bike_profile";

        List<Profile> profiles = asList(
                new Profile(bikeProfile).setVehicle("bike")
        );

        GraphHopper hopperCustomSpeeds = new GraphHopperCustomSpeeds(new SpeedsTooHighSpeedProvider()).
                setGraphHopperLocation(GH_LOCATION_CUSTOM_SPEEDS_RELOAD).
                setOSMFile(MONACO).
                setProfiles(profiles).
                setStoreOnFlush(true)
                .setEncodedValuesString("roundabout, road_class, road_class_link, road_environment, max_speed, road_access, ferry_speed, bike_network, get_off_bike, smoothness, osm_way_id");
        hopperCustomSpeeds.getCHPreparationHandler().setCHProfiles(
                new CHProfile(bikeProfile)
        );

        hopperCustomSpeeds.importOrLoad();

        double maxSpeed = hopperCustomSpeeds.getEncodingManager().getDecimalEncodedValue(VehicleSpeed.key("bike")).getMaxOrMaxStorableDecimal();

       hopperCustomSpeeds.close();

        GraphHopper loadHopperCustomSpeeds = new GraphHopper().
                setGraphHopperLocation(GH_LOCATION_CUSTOM_SPEEDS_RELOAD).
                setOSMFile(MONACO).
                setProfiles(profiles)
                .setEncodedValuesString("roundabout, road_class, road_class_link, road_environment, max_speed, road_access, ferry_speed, bike_network, get_off_bike, smoothness, osm_way_id");
        hopperCustomSpeeds.getCHPreparationHandler().setCHProfiles(
                new CHProfile(bikeProfile)
        );

        loadHopperCustomSpeeds.importOrLoad();

        double loadMaxSpeed = loadHopperCustomSpeeds.getEncodingManager().getDecimalEncodedValue(VehicleSpeed.key("bike")).getMaxOrMaxStorableDecimal();

        loadHopperCustomSpeeds.close();

        assertEquals(maxSpeed, loadMaxSpeed, 1);

    }



}

