package com.graphhopper.speeds;

import java.util.Optional;

public interface WaySpeedsProvider {
    double speedForWay(long osmWayId);
}
