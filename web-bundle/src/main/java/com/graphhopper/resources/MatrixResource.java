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
package com.graphhopper.resources;

import com.fasterxml.jackson.databind.node.JsonNodeFactory;
import com.fasterxml.jackson.databind.node.ObjectNode;
import com.graphhopper.routing.matrix.GHMatrixRequest;
import com.graphhopper.routing.matrix.GHMatrixResponse;
import com.graphhopper.routing.matrix.DistanceMatrix;
import com.graphhopper.GraphHopper;
import com.graphhopper.http.GHPointParam;
import com.graphhopper.http.GHRequestTransformer;
import com.graphhopper.http.ProfileResolver;
import com.graphhopper.jackson.MultiException;
import com.graphhopper.jackson.ResponsePathSerializer;
import com.graphhopper.util.*;
import com.graphhopper.util.shapes.GHPoint;
import io.dropwizard.jersey.params.AbstractParam;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.inject.Inject;
import javax.inject.Named;
import javax.servlet.http.HttpServletRequest;
import javax.validation.constraints.NotNull;
import javax.ws.rs.*;
import javax.ws.rs.core.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import static com.graphhopper.util.Parameters.Details.PATH_DETAILS;
import static com.graphhopper.util.Parameters.Routing.*;
import static java.util.stream.Collectors.toList;

/**
 * Resource to use GraphHopper in a remote client application like mobile or browser. Note: If type
 * is json it returns the points in GeoJson array format [longitude,latitude] unlike the format "lat,lon"
 * used for the request. See the full API response format in docs/web/api-doc.md
 *
 * @author Peter Karich
 */
@Path("matrix")
public class MatrixResource {

    private static final Logger logger = LoggerFactory.getLogger(MatrixResource.class);

    private final GraphHopper graphHopper;
    private final ProfileResolver profileResolver;
    private final GHRequestTransformer ghRequestTransformer;
    private final Boolean hasElevation;
    private final String osmDate;

    @Inject
    public MatrixResource(GraphHopper graphHopper, ProfileResolver profileResolver, GHRequestTransformer ghRequestTransformer, @Named("hasElevation") Boolean hasElevation) {
        this.graphHopper = graphHopper;
        this.profileResolver = profileResolver;
        this.ghRequestTransformer = ghRequestTransformer;
        this.hasElevation = hasElevation;
        this.osmDate = graphHopper.getProperties().get("datareader.data.date");
    }

    @GET
    @Produces(MediaType.APPLICATION_JSON)
    public Response doGet(
            @Context HttpServletRequest httpReq,
            @Context UriInfo uriInfo,
            @QueryParam("type") @DefaultValue("json") String type,
            @QueryParam(WAY_POINT_MAX_DISTANCE) @DefaultValue("1") double minPathPrecision,
            @QueryParam(ALGORITHM) @DefaultValue("") String algoStr,
            @QueryParam("profile") String profileName,
            @QueryParam("point") @NotNull List<GHPointParam> pointParams,
            @QueryParam("from_point") @NotNull List<GHPointParam> fromPointParams,
            @QueryParam("to_point") @NotNull List<GHPointParam> toPointParams,
            @QueryParam(POINT_HINT) List<String> pointHints,
            @QueryParam("from_point_hint") List<String> fromPointHints,
            @QueryParam("to_point_hint") List<String> toPointHints,
            @QueryParam(SNAP_PREVENTION) List<String> snapPreventions,
            @QueryParam(CURBSIDE) List<String> curbsides,
            @QueryParam("from_curbside") List<String> fromCurbsides,
            @QueryParam("to_curbside") List<String> toCurbsides,
            @QueryParam("out_array") List<String> outArrays,
            @QueryParam("fail_fast") @DefaultValue("true") String failFast) {
        StopWatch sw = new StopWatch().start();
        List<GHPoint> points = pointParams.stream().map(AbstractParam::get).collect(toList());
        List<GHPoint> fromPoints = fromPointParams.stream().map(AbstractParam::get).collect(toList());
        List<GHPoint> toPoints = toPointParams.stream().map(AbstractParam::get).collect(toList());

        GHMatrixRequest request = new GHMatrixRequest();
        initHints(request.getHints(), uriInfo.getQueryParameters());

        request.setAlgorithm(algoStr).
                setProfile(profileName);
        if (points != null) {
            request.setOrigins(points).
                    setDestinations(points);
        } else {
            request.setOrigins(fromPoints).
                    setDestinations(toPoints);
        }
        if (pointHints != null) {
            request.setOriginPointHints(pointHints).
                    setDestinationPointHints(pointHints);
        } else {
            request.setOriginPointHints(fromPointHints).
                    setDestinationPointHints(toPointHints);
        }
        request.setSnapPreventions(snapPreventions);
        if (curbsides != null) {
            request.setOriginCurbsides(curbsides).
                    setDestinationCurbsides(curbsides);
        } else {
            request.setOriginCurbsides(fromCurbsides).
                    setDestinationCurbsides(toCurbsides);
        }
        request.setFailFast(Boolean.parseBoolean(failFast));

        // request = ghRequestTransformer.transformRequest(request);

        PMap profileResolverHints = new PMap(request.getHints());
        profileResolverHints.putObject("profile", profileName);
        profileName = profileResolver.resolveProfile(profileResolverHints);
        removeLegacyParameters(request.getHints());
        request.setProfile(profileName);

        GHMatrixResponse ghResponse = graphHopper.matrix(request);

        double took = sw.stop().getMillisDouble();
        String logStr = (httpReq.getRemoteAddr() + " " + httpReq.getLocale() + " " + httpReq.getHeader("User-Agent")) + " " + fromPoints + " " + toPoints + ", took: " + String.format("%.1f", took) + "ms, algo: " + algoStr + ", profile: " + profileName;
        logger.info(logStr);

        DistanceMatrix matrix = ghResponse.getMatrix();
        ObjectNode json = JsonNodeFactory.instance.objectNode();
        final ObjectNode info = json.putObject("info");
        info.putPOJO("copyrights", ResponsePathSerializer.COPYRIGHTS);
        info.put("took", Math.round((float) sw.getMillis()));
        if (outArrays != null && outArrays.contains("times")) {
            json.putPOJO("times", matrix.getTimes());
        }
        if (outArrays != null && outArrays.contains("distances")) {
            json.putPOJO("distances", matrix.getDistances());
        }

        return Response.ok(json).
                header("X-GH-Took", "" + Math.round(took)).
                type(MediaType.APPLICATION_JSON).
                build();
    }

    @POST
    @Consumes(MediaType.APPLICATION_JSON)
    @Produces(MediaType.APPLICATION_JSON)
    public Response doPost(@NotNull GHMatrixRequest request, @Context HttpServletRequest httpReq) {
        StopWatch sw = new StopWatch().start();
        // request = ghRequestTransformer.transformRequest(request);

        PMap profileResolverHints = new PMap(request.getHints());
        profileResolverHints.putObject("profile", request.getProfile());
        request.setProfile(profileResolver.resolveProfile(profileResolverHints));
        removeLegacyParameters(request.getHints());

        GHMatrixResponse ghResponse = graphHopper.matrix(request);
        List<String> outArrays = request.getHints().getObject("out_arrays", new ArrayList<String>());

        double took = sw.stop().getMillisDouble();
        String infoStr = httpReq.getRemoteAddr() + " " + httpReq.getLocale() + " " + httpReq.getHeader("User-Agent");
        String logStr = infoStr + " " + request.getOrigins().size() + "x" + request.getDestinations().size() + ", took: "
                + String.format("%.1f", took) + " ms, algo: " + request.getAlgorithm() + ", profile: " + request.getProfile()
                + ", custom_model: " + request.getCustomModel();
        logger.info(logStr);

        DistanceMatrix matrix = ghResponse.getMatrix();
        ObjectNode json = JsonNodeFactory.instance.objectNode();
        final ObjectNode info = json.putObject("info");
        info.putPOJO("copyrights", ResponsePathSerializer.COPYRIGHTS);
        info.put("took", Math.round((float) sw.getMillis()));
        if (outArrays != null && outArrays.contains("times")) {
            json.putPOJO("times", matrix.getTimes());
        }
        if (outArrays != null && outArrays.contains("distances")) {
            json.putPOJO("distances", matrix.getDistances());
        }

        return Response.ok(json).
                header("X-GH-Took", "" + Math.round(took)).
                type(MediaType.APPLICATION_JSON).
                build();
    }

    public static void removeLegacyParameters(PMap hints) {
        // these parameters should only be used to resolve the profile, but should not be passed to GraphHopper
        hints.remove("weighting");
        hints.remove("vehicle");
        hints.remove("edge_based");
        hints.remove("turn_costs");
    }

    static void initHints(PMap m, MultivaluedMap<String, String> parameterMap) {
        for (Map.Entry<String, List<String>> e : parameterMap.entrySet()) {
            if (e.getValue().size() == 1) {
                m.putObject(Helper.camelCaseToUnderScore(e.getKey()), Helper.toObject(e.getValue().get(0)));
            } else {
                // TODO e.g. 'point' parameter occurs multiple times and we cannot throw an exception here
                //  unknown parameters (hints) should be allowed to be multiparameters, too, or we shouldn't use them for
                //  known parameters either, _or_ known parameters must be filtered before they come to this code point,
                //  _or_ we stop passing unknown parameters altogether.
                // throw new WebApplicationException(String.format("This query parameter (hint) is not allowed to occur multiple times: %s", e.getKey()));
                // see also #1976
            }
        }
    }
}
