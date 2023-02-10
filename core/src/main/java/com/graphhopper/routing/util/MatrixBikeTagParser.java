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

import com.graphhopper.reader.ReaderWay;
import com.graphhopper.routing.ev.*;
import com.graphhopper.storage.IntsRef;
import com.graphhopper.util.PMap;

import static com.graphhopper.routing.ev.RouteNetwork.*;
import static com.graphhopper.routing.ev.RouteNetwork.LOCAL;
import static com.graphhopper.routing.util.PriorityCode.*;

/**
 * Specifies the settings for cycletouring/trekking
 *
 * @author ratrun
 * @author Peter Karich
 */
public class MatrixBikeTagParser extends BikeCommonTagParser {

    protected static final int PUSHING_SECTION_SPEED = 5;

    public MatrixBikeTagParser(EncodedValueLookup lookup, PMap properties) {
        this(
                lookup.getBooleanEncodedValue(VehicleAccess.key(properties.getString("name", "matrixbike"))),
                lookup.getDecimalEncodedValue(VehicleSpeed.key(properties.getString("name", "matrixbike"))),
                lookup.getDecimalEncodedValue(VehiclePriority.key(properties.getString("name", "matrixbike"))),
                lookup.getEnumEncodedValue(BikeNetwork.KEY, RouteNetwork.class),
                lookup.getEnumEncodedValue(Smoothness.KEY, Smoothness.class),
                properties.getString("name", "matrixbike"),
                lookup.getBooleanEncodedValue(Roundabout.KEY),
                lookup.hasEncodedValue(TurnCost.key(properties.getString("name", "matrixbike"))) ? lookup.getDecimalEncodedValue(TurnCost.key(properties.getString("name", "matrixbike"))) : null
        );
        blockPrivate(properties.getBool("block_private", true));
        blockFords(properties.getBool("block_fords", false));
    }

    public MatrixBikeTagParser(BooleanEncodedValue accessEnc, DecimalEncodedValue speedEnc, DecimalEncodedValue priorityEnc,
                               EnumEncodedValue<RouteNetwork> bikeRouteEnc, EnumEncodedValue<Smoothness> smoothnessEnc,
                               String name, BooleanEncodedValue roundaboutEnc, DecimalEncodedValue turnCostEnc) {
        super(accessEnc, speedEnc, priorityEnc, bikeRouteEnc, smoothnessEnc, name, roundaboutEnc, turnCostEnc);
        addPushingSection("path");
        addPushingSection("footway");
        addPushingSection("pedestrian");
        addPushingSection("steps");
        addPushingSection("platform");
        addPushingSection("pier");

        /* START */

        setTrackTypeSpeed("grade1", 11); // paved
        setTrackTypeSpeed("grade2", 11); // now unpaved ...
        setTrackTypeSpeed("grade3", 8);
        setTrackTypeSpeed("grade4", 6);
        setTrackTypeSpeed("grade5", 4); // like sand/grass

        setSurfaceSpeed("paved", 11);
        setSurfaceSpeed("asphalt", 11);
        setSurfaceSpeed("cobblestone", 11);
        setSurfaceSpeed("cobblestone:flattened", 11);
        setSurfaceSpeed("sett", 11);
        setSurfaceSpeed("concrete", 11);
        setSurfaceSpeed("concrete:lanes", 11);
        setSurfaceSpeed("concrete:plates", 11);
        setSurfaceSpeed("paving_stones", 11);
        setSurfaceSpeed("paving_stones:30", 11);
        setSurfaceSpeed("unpaved", 11);
        setSurfaceSpeed("compacted", 11);
        setSurfaceSpeed("dirt", 11);
        setSurfaceSpeed("earth", 11);
        setSurfaceSpeed("fine_gravel", 11);
        setSurfaceSpeed("grass", 11);
        setSurfaceSpeed("grass_paver", 11);
        setSurfaceSpeed("gravel", 11);
        setSurfaceSpeed("ground", 11);
        setSurfaceSpeed("metal", 11);
        setSurfaceSpeed("mud", 11);
        setSurfaceSpeed("pebblestone", 11);
        setSurfaceSpeed("salt", 11);
        setSurfaceSpeed("sand", 11);
        setSurfaceSpeed("wood", 11);

        setHighwaySpeed("living_street", 11);
        setHighwaySpeed("steps", 2);
        avoidHighwayTags.add("steps");

        setHighwaySpeed("cycleway", 11);
        setHighwaySpeed("path", 10);
        setHighwaySpeed("footway", 11);
        setHighwaySpeed("platform", PUSHING_SECTION_SPEED);
        setHighwaySpeed("pedestrian", 11);
        setHighwaySpeed("track", 11);
        setHighwaySpeed("service", 11);
        setHighwaySpeed("residential", 11);
        // no other highway applies:
        setHighwaySpeed("unclassified", 11);
        // unknown road:
        setHighwaySpeed("road", 11);

        setHighwaySpeed("trunk", 11);
        setHighwaySpeed("trunk_link", 11);
        setHighwaySpeed("primary", 11);
        setHighwaySpeed("primary_link", 11);
        setHighwaySpeed("secondary", 11);
        setHighwaySpeed("secondary_link", 11);
        setHighwaySpeed("tertiary", 11);
        setHighwaySpeed("tertiary_link", 11);

        // special case see tests and #191
        setHighwaySpeed("motorway", 11);
        setHighwaySpeed("motorway_link", 11);
        avoidHighwayTags.add("motorway");
        avoidHighwayTags.add("motorway_link");

        setHighwaySpeed("bridleway", PUSHING_SECTION_SPEED);
        avoidHighwayTags.add("bridleway");
        avoidHighwayTags.add("trunk");
        avoidHighwayTags.add("trunk_link");
        avoidHighwayTags.add("primary");
        avoidHighwayTags.add("primary_link");
        avoidHighwayTags.add("secondary");
        avoidHighwayTags.add("secondary_link");

        // preferHighwayTags.add("road");
        preferHighwayTags.add("service");
        preferHighwayTags.add("tertiary");
        preferHighwayTags.add("tertiary_link");
        preferHighwayTags.add("residential");
        preferHighwayTags.add("unclassified");

        barriers.add("kissing_gate");
        barriers.add("stile");
        barriers.add("turnstile");

        setSpecificClassBicycle("touring");
    }

    @Override
    public IntsRef handleWayTags(IntsRef edgeFlags, ReaderWay way) {
        WayAccess access = getAccess(way);

        if (access.canSkip())
            return edgeFlags;

        Integer priorityFromRelation = routeMap.get(bikeRouteEnc.getEnum(false, edgeFlags));
        double wayTypeSpeed = getSpeed(way);
        if (!access.isFerry()) {
            wayTypeSpeed = applyMaxSpeed(way, wayTypeSpeed);
            Smoothness smoothness = smoothnessEnc.getEnum(false, edgeFlags);
            wayTypeSpeed = Math.max(MIN_SPEED, smoothnessFactor.get(smoothness) * wayTypeSpeed);

            avgSpeedEnc.setDecimal(false, edgeFlags, wayTypeSpeed);
            if (avgSpeedEnc.isStoreTwoDirections())
                avgSpeedEnc.setDecimal(true, edgeFlags, wayTypeSpeed);
            handleAccess(edgeFlags, way);
            handlerBikePush(edgeFlags, way);
        } else {
            double ferrySpeed = ferrySpeedCalc.getSpeed(way);
            avgSpeedEnc.setDecimal(false, edgeFlags, ferrySpeed);
            if (avgSpeedEnc.isStoreTwoDirections())
                avgSpeedEnc.setDecimal(true, edgeFlags, ferrySpeed);
            accessEnc.setBool(false, edgeFlags, true);
            accessEnc.setBool(true, edgeFlags, true);
            priorityFromRelation = SLIGHT_AVOID.getValue();
        }

        priorityEnc.setDecimal(false, edgeFlags, PriorityCode.getValue(handlePriority(way, wayTypeSpeed, priorityFromRelation)));
        return edgeFlags;
    }

    private void handlerBikePush(IntsRef edgeFlags,ReaderWay way){

        boolean backwardInaccessible = !accessEnc.getBool(true,edgeFlags);
        boolean forwardInaccessible = !accessEnc.getBool(false,edgeFlags);

        double backwardSpeed = avgSpeedEnc.getDecimal(true,edgeFlags);
        double forwardSpeed = avgSpeedEnc.getDecimal(false,edgeFlags);

        String highwayTag = way.getTag("highway");
        double pushHighwaySpeed = PUSHING_SECTION_SPEED;
        if(highwaySpeeds.containsKey(highwayTag)){
            double highwaySpeed = highwaySpeeds.get(highwayTag);
            if(highwaySpeed < pushHighwaySpeed) {
                pushHighwaySpeed = highwaySpeed;
            }
        }

        // pushing bikes - if no other mode found
        if(forwardInaccessible || backwardInaccessible
                || forwardSpeed == Double.POSITIVE_INFINITY
                || backwardSpeed == Double.POSITIVE_INFINITY){

            if(!way.hasTag("foot", "no")){

                boolean implyOneWay = way.hasTag("junction", "roundabout")
                        || way.hasTag("junction", "circular")
                        || way.hasTag("highway", "motorway");

                boolean wayTypeAllowPushing = way.hasTag("railway","platform")
                        || way.hasTag("bridge","movable")
                        || way.hasTag("public_transport","platform")
                        || way.hasTag("amenity","parking","parking_entrance")
                        || highwaySpeeds.containsKey(highwayTag)
                        || way.hasTag("access","yes","permissive","designated"); //intendedValues?

                double pushForwardSpeed = Double.POSITIVE_INFINITY;
                double pushBackwardSpeed = Double.POSITIVE_INFINITY;

                if (way.hasTag("highway", pushingSectionsHighways)){
                    pushForwardSpeed = pushHighwaySpeed;
                    pushBackwardSpeed = pushHighwaySpeed;
                }else{
                    if(way.hasTag("foot", "yes")){

                        pushForwardSpeed = pushHighwaySpeed;
                        if(!implyOneWay){
                            pushBackwardSpeed = pushHighwaySpeed;
                        }
                    }else if(way.hasTag("foot:forward", "yes")){
                        pushForwardSpeed = pushHighwaySpeed;
                    }else if(way.hasTag("foot:backward", "yes")) {
                        pushBackwardSpeed = pushHighwaySpeed;
                    }else if(wayTypeAllowPushing) {
                        pushForwardSpeed = pushHighwaySpeed;
                        if(!implyOneWay){
                            pushBackwardSpeed = pushHighwaySpeed;
                        }
                    }
                }

                if(pushForwardSpeed != Double.POSITIVE_INFINITY && (forwardInaccessible || forwardSpeed == Double.POSITIVE_INFINITY)){
                    accessEnc.setBool(false, edgeFlags, true);
                    avgSpeedEnc.setDecimal(false, edgeFlags, pushForwardSpeed);
                }

                if(pushBackwardSpeed != Double.POSITIVE_INFINITY && (backwardInaccessible || backwardSpeed == Double.POSITIVE_INFINITY)){
                    if(accessEnc.isStoreTwoDirections() && avgSpeedEnc.isStoreTwoDirections()) {
                        accessEnc.setBool(true, edgeFlags, true);
                        avgSpeedEnc.setDecimal(true, edgeFlags, pushBackwardSpeed);
                    }
                }
            }
        }

        // dismount
        if (way.hasTag("bicycle", "dismount")){
            accessEnc.setBool(false, edgeFlags, true);
            avgSpeedEnc.setDecimal(false, edgeFlags, PUSHING_SECTION_SPEED);

            if(accessEnc.isStoreTwoDirections() && avgSpeedEnc.isStoreTwoDirections()){
                accessEnc.setBool(true, edgeFlags, true);
                avgSpeedEnc.setDecimal(true, edgeFlags, PUSHING_SECTION_SPEED);
            }
        }

    }

}
