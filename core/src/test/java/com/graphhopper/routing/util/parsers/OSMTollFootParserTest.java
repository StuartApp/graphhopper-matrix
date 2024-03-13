package com.graphhopper.routing.util.parsers;

import com.graphhopper.reader.ReaderWay;
import com.graphhopper.routing.ev.EncodedValue;
import com.graphhopper.routing.ev.EnumEncodedValue;
import com.graphhopper.routing.ev.TollFoot;
import com.graphhopper.storage.IntsRef;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class OSMTollFootParserTest {
    private EnumEncodedValue<TollFoot> tollEnc;
    private OSMTollFootParser parser;

    @BeforeEach
    public void setUp() {
        tollEnc = new EnumEncodedValue<>(TollFoot.KEY, TollFoot.class);
        tollEnc.init(new EncodedValue.InitializerConfig());
        parser = new OSMTollFootParser(tollEnc);
    }

    @Test
    public void testSimpleTags() {
        ReaderWay readerWay = new ReaderWay(1);
        IntsRef relFlags = new IntsRef(2);
        IntsRef intsRef = new IntsRef(1);
        readerWay.setTag("highway", "primary");
        parser.handleWayTags(intsRef, readerWay, relFlags);
        assertEquals(TollFoot.MISSING, tollEnc.getEnum(false, intsRef));

        intsRef = new IntsRef(1);
        readerWay.setTag("highway", "primary");
        readerWay.setTag("toll:foot", "yes");
        parser.handleWayTags(intsRef, readerWay, relFlags);
        assertEquals(TollFoot.YES, tollEnc.getEnum(false, intsRef));

        intsRef = new IntsRef(1);
        readerWay.setTag("highway", "primary");
        readerWay.setTag("toll:foot", "no");
        parser.handleWayTags(intsRef, readerWay, relFlags);
        assertEquals(TollFoot.NO, tollEnc.getEnum(false, intsRef));
    }
}