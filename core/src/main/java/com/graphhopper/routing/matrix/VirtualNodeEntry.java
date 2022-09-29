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

public class VirtualNodeEntry{

    public int node;
    public int edge;
    public double weight;
    public double distance;
    public long time;
    public boolean virtual;
    public boolean combined;
    public int traversalId;


    public VirtualNodeEntry(int traversalId, int node, int edge, double weight, long time, double distance) {
        this.traversalId = traversalId;
        this.node = node;
        this.edge = edge;
        this.weight = weight;
        this.time = time;
        this.distance = distance;
        this.combined = false;
    }

    public VirtualNodeEntry(int traversalId,int node, int edge, double weight, long time, double distance,boolean combined) {
        this(traversalId,node,edge,weight,time,distance);
        this.combined = combined;
    }

    public VirtualNodeEntry combine(VirtualNodeEntry other) {
        return new VirtualNodeEntry(this.traversalId,this.node,this.edge,this.weight + other.weight, this.time + other.time,
                this.distance + other.distance);
    }

    public double getWeightOfVisitedPath() {
        return weight;
    }

    @Override
    public String toString() {
        return node + " weight: " + weight + " time: " + time + " distance :" + distance;
    }

}
