package com.graphhopper.speeds;

import com.graphhopper.routing.ev.RoadClass;

import java.util.Optional;

public interface WaySpeedsProvider {
    double speedForWay(long osmWayId, RoadClass roadClass);
}
