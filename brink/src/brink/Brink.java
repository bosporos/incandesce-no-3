package brink;

import processing.core.*;
import processing.data.*;
import processing.event.*;
import processing.opengl.*;

/**
 * Basically just a container for a PApplet instance.
 */
public class Brink {

    /**
     * The PApplet instance.
     */
    public static PApplet applet_instance;

    /**
     * We need an instance of PApplet for a lot of Processing's functionality
     */
    public static void init(PApplet parent) {
        Brink.applet_instance = parent;
    }
}
