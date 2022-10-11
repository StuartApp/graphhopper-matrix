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
package com.graphhopper.routing.matrix;

import com.graphhopper.util.EdgeIterator;
import com.graphhopper.util.EdgeIteratorState;

import java.util.Objects;


public class MatrixEntry implements Comparable<MatrixEntry> {

    public int edge;
    public int adjNode;
    public int baseNode;
    public double weight;
    public double distance;
    public long time;
    public int level;

    /**
     * The id of the incoming original edge at this shortest path tree entry. For original edges this is the same
     * as the edge id, but for shortcuts this is the id of the last original edge of the shortcut.
     *
     * @see EdgeIteratorState#getOrigEdgeLast()
     */
    public int originalEdge;


    public MatrixEntry(int node, double weight, long time, double distance, int level) {
        this(EdgeIterator.NO_EDGE,EdgeIterator.NO_EDGE, node,node, level, weight,time,distance);
    }

    public MatrixEntry(int edge, int origEdgeId,int adjNode, int baseNode, int level, double weight, long time, double distance) {
        this.edge = edge;
        this.originalEdge = origEdgeId;
        this.adjNode = adjNode;
        this.baseNode = baseNode;
        this.weight = weight;
        this.time = time;
        this.distance = distance;
        this.level = level;
    }

    public double getWeightOfVisitedPath() {
        return weight;
    }

    @Override
    public String toString() {
        return "MatrixEntry{" +
                "edge=" + edge +
                ", adjNode=" + adjNode +
                ", baseNode=" + baseNode +
                ", weight=" + weight +
                ", distance=" + distance +
                ", time=" + time +
                ", level=" + level +
                ", originalEdge=" + originalEdge +
                '}';
    }

    @Override
    public int compareTo(MatrixEntry o) {

        //if(level < o.level){
            //return -1;
        //}else {
            if (level < o.level)
                return -1;

            // assumption no NaN and no -0
            return level > o.level ? 1 : 0;
        //}

    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        MatrixEntry that = (MatrixEntry) o;
        return edge == that.edge && adjNode == that.adjNode && baseNode == that.baseNode && originalEdge == that.originalEdge;
    }

    @Override
    public int hashCode() {
        return Objects.hash(edge, adjNode, baseNode, originalEdge);
    }
}
