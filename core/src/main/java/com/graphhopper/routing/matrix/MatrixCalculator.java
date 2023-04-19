package com.graphhopper.routing.matrix;

import com.graphhopper.storage.index.Snap;

import java.util.List;

public interface MatrixCalculator {
    DistanceMatrix calcMatrix(MatrixSnapResult origins, MatrixSnapResult destinations);

    String getDebugString();

    int getVisitedNodes();
}
