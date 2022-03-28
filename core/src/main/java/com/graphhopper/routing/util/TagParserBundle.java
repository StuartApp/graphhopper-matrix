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

package com.graphhopper.routing.util;

import com.graphhopper.reader.ReaderRelation;
import com.graphhopper.reader.ReaderWay;
import com.graphhopper.routing.ev.EncodedValue;
import com.graphhopper.routing.util.parsers.RelationTagParser;
import com.graphhopper.routing.util.parsers.TagParser;
import com.graphhopper.routing.util.parsers.TurnCostParser;
import com.graphhopper.storage.IntsRef;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

public class TagParserBundle {
    private final List<TagParser> wayTagParsers;
    private final List<RelationTagParser> relationTagParsers;
    private final List<TurnCostParser> turnCostParsers;
    private final EncodedValue.InitializerConfig relConfig = new EncodedValue.InitializerConfig();

    public TagParserBundle() {
        this(new ArrayList<>(), new ArrayList<>(), new ArrayList<>());
    }

    public TagParserBundle(List<TagParser> wayTagParsers, List<RelationTagParser> relationTagParsers, List<TurnCostParser> turnCostParsers) {
        this.wayTagParsers = wayTagParsers;
        this.relationTagParsers = relationTagParsers;
        this.turnCostParsers = turnCostParsers;
    }

    public TagParserBundle addWayTagParser(TagParser tagParser) {
        wayTagParsers.add(tagParser);
        return this;
    }

    public TagParserBundle addRelationTagParser(Function<EncodedValue.InitializerConfig, RelationTagParser> createRelationTagParser) {
        relationTagParsers.add(createRelationTagParser.apply(relConfig));
        return this;
    }

    public TagParserBundle addTurnCostTagParser(TurnCostParser turnCostParser) {
        turnCostParsers.add(turnCostParser);
        return this;
    }

    public IntsRef handleRelationTags(ReaderRelation relation, IntsRef relFlags) {
        for (RelationTagParser relParser : relationTagParsers) {
            relParser.handleRelationTags(relFlags, relation);
        }
        return relFlags;
    }

    public IntsRef handleWayTags(IntsRef edgeFlags, ReaderWay way, IntsRef relationFlags) {
        for (RelationTagParser relParser : relationTagParsers)
            relParser.handleWayTags(edgeFlags, way, relationFlags);
        for (TagParser parser : wayTagParsers)
            parser.handleWayTags(edgeFlags, way, relationFlags);
        return edgeFlags;
    }

    public IntsRef createRelationFlags() {
        // todonow: we can determine this from relConfig
        return new IntsRef(2);
    }
}
