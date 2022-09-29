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


/*
 * Model: (Edge) ---->  (Adjacency Node) -----> (Next Edge)
 */

public class RankedNode implements Comparable<RankedNode> {
    public int edge;
    public int adjNode;
    public int rank;
    public int traversalId;
    public boolean virtual;


    public RankedNode(int traversalId,int edge,int adjNode,int rank, boolean virtual) {
        this.traversalId = traversalId;
        this.edge = edge;
        this.adjNode = adjNode;
        this.rank = rank;
        this.virtual = virtual;
    }


    @Override
    public int compareTo(RankedNode o) {

        if (rank < o.rank)
            return -1;

        // assumption no NaN and no -0
        return rank > o.rank ? 1 : 0;
    }

    @Override
    public String toString() {
        return "--" + edge + "--> " + adjNode  + " rank: " + rank + " traversal: " + traversalId;
    }
}
