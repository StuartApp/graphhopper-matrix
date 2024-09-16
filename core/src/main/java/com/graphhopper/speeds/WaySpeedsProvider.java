package com.graphhopper.speeds;

import com.graphhopper.routing.ev.RoadClass;

import java.util.Optional;

public interface WaySpeedsProvider {
    Optional<SpeedKmByHour> speedForWay(long osmWayId);
    Optional<SpeedKmByHour> speedForRoadClass(RoadClass roadClass);
    void setMetrics(Metrics metrics);
    Metrics getMetrics();
}
