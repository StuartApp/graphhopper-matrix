package com.graphhopper.speeds;

import com.graphhopper.routing.ev.RoadClass;

import java.util.Optional;

public interface WaySpeedsProvider {
    Optional<Double> speedForWay(long osmWayId);
    Optional<Double> speedForRoadClass(RoadClass roadClass);
}
