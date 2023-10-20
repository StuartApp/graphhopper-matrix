package com.graphhopper.routing.util.parsers;

import com.graphhopper.routing.ev.*;
import com.graphhopper.util.PMap;

public class MatrixBikeAverageSpeedParser extends BikeCommonAverageSpeedParser {

    protected static final int PUSHING_SECTION_SPEED = 5;
    protected static final int MIN_SPEED = 2;

    public MatrixBikeAverageSpeedParser(EncodedValueLookup lookup, PMap properties) {
        this(lookup.getDecimalEncodedValue(VehicleSpeed.key(properties.getString("name", "matrixbike"))),
                lookup.getEnumEncodedValue(Smoothness.KEY, Smoothness.class),lookup.getDecimalEncodedValue(FerrySpeed.KEY));
    }

    public MatrixBikeAverageSpeedParser(DecimalEncodedValue speedEnc, EnumEncodedValue<Smoothness> smoothnessEnc, DecimalEncodedValue ferrySpeedEnc) {
        super(speedEnc, smoothnessEnc, ferrySpeedEnc);


        //Custom
        addPushingSection("path");
        addPushingSection("footway");
        addPushingSection("pedestrian");
        addPushingSection("steps");
        addPushingSection("platform");
        addPushingSection("pier");

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
        setSurfaceSpeed("ice", MIN_SPEED);
        setSurfaceSpeed("metal", 11);
        setSurfaceSpeed("mud", 11);
        setSurfaceSpeed("pebblestone", 11);
        setSurfaceSpeed("salt", PUSHING_SECTION_SPEED);
        setSurfaceSpeed("sand", PUSHING_SECTION_SPEED);
        setSurfaceSpeed("wood", PUSHING_SECTION_SPEED);

        setHighwaySpeed("living_street", 11);
        setHighwaySpeed("steps", 2);

        setHighwaySpeed("cycleway", 11);
        setHighwaySpeed("path", 10);
        setHighwaySpeed("footway", 11);
        setHighwaySpeed("platform", PUSHING_SECTION_SPEED);
        setHighwaySpeed("pedestrian", PUSHING_SECTION_SPEED);
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

        setHighwaySpeed("bridleway", PUSHING_SECTION_SPEED);

        // note that this factor reduces the speed but only until MIN_SPEED
        setSmoothnessSpeedFactor(Smoothness.MISSING, 1.0d);
        setSmoothnessSpeedFactor(Smoothness.OTHER, 0.7d);
        setSmoothnessSpeedFactor(Smoothness.EXCELLENT, 1.1d);
        setSmoothnessSpeedFactor(Smoothness.GOOD, 1.0d);
        setSmoothnessSpeedFactor(Smoothness.INTERMEDIATE, 0.9d);
        setSmoothnessSpeedFactor(Smoothness.BAD, 0.7d);
        setSmoothnessSpeedFactor(Smoothness.VERY_BAD, 0.4d);
        setSmoothnessSpeedFactor(Smoothness.HORRIBLE, 0.3d);
        setSmoothnessSpeedFactor(Smoothness.VERY_HORRIBLE, 0.1d);
        setSmoothnessSpeedFactor(Smoothness.IMPASSABLE, 0);

        intendedValues.add("yes");
        intendedValues.add("designated");
        intendedValues.add("official");
        intendedValues.add("permissive");
        //End Custom
    }
}
