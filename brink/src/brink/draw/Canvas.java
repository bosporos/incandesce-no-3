package brink.draw;

import brink.Brink;

import processing.core.PApplet;
import processing.core.PConstants;

public class Canvas implements PConstants {
    private static Canvas instance = null;
    public static Canvas instance() {
        if(Canvas.instance == null) {
            Canvas.instance = new Canvas();
        }
        return Canvas.instance;
    }

    public PApplet inner;

    private Canvas() {
        inner = Brink.applet_instance;
    }

    public void init() {
        this.inner.colorMode(HSB, 356, 100, 100, 1);
        this.inner.ellipseMode(RADIUS);
    }

    public void fill(Color color) {
        this.inner.fill(color.inner);
    }

    public void stroke(Color color) {
        this.inner.stroke(color.inner);
    }

    public void background(Color color) {
        this.inner.background(color.inner);
    }
}
