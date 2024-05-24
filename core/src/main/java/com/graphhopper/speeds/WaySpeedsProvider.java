package com.graphhopper.speeds;

import com.graphhopper.routing.ev.RoadClass;

import java.util.Optional;

public interface WaySpeedsProvider {
    Optional<Double> speedKmHourForWay(long osmWayId);
    Optional<Double> speedKmHourForRoadClass(RoadClass roadClass);
}
