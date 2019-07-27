package brink.draw;

import brink.Brink;

public class Color {
    public int inner;

    public Color(float h, float s, float b) {
        this.inner = Brink.applet_instance.color(h, s, b, 1);
    }

    public Color(float h, float s, float b, float a) {
        this.inner = Brink.applet_instance.color(h, s, b, a);
    }

    public Color(Color origin) {
        this.inner = origin.inner;
    }

    public Color(int origin) {
        this.inner = origin;
    }

    public Color set_alpha(float a) {
        this.inner = ((int)(a * 0xFF) << 24) | (this.inner & 0x00FFFFFF);
        return this;
    }

    @Override
    public Color clone() {
        return new Color(this);
    }

    public Color apply_alpha(float a) {
        return this.clone().set_alpha(a);
    }

}
