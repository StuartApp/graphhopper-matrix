package com.graphhopper.routing.matrix;

import com.carrotsearch.hppc.IntArrayList;
import com.graphhopper.storage.index.Snap;

import java.util.List;


public class MatrixSnapResult {

    List<Snap> snaps;
    IntArrayList pointsNotFound = new IntArrayList();

    public MatrixSnapResult(List<Snap> snaps, IntArrayList pointsNotFound) {
        this.snaps = snaps;
        this.pointsNotFound = pointsNotFound;
    }
}
