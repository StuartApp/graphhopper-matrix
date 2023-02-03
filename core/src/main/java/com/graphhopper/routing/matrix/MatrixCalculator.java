package com.graphhopper.routing.matrix;

import com.graphhopper.storage.index.Snap;

import java.util.List;

public interface MatrixCalculator {
    DistanceMatrix calcMatrix(List<Snap> origins, List<Snap> destinations);
    DistanceMatrix calcMatrixV2(MatrixSnapResult origins, MatrixSnapResult destinations);

    String getDebugString();

    int getVisitedNodes();
}
