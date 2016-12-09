package vis;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.geom.Point2D;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

import no.geosoft.cc.geometry.Geometry;
import no.geosoft.cc.geometry.Matrix4x4;
import no.geosoft.cc.graphics.*;

public class Vis extends JFrame
        implements ActionListener, GInteraction{

    private JButton  frontButton_;
    private JButton  backButton_;
    private JButton  forwardButton_;
    private JButton  backwardButton_;

    private GScene   scene_;
    private GObject  interactionObject_;
    private Color    color_;
    private boolean drawFlag = true;
    private int timeStep;
    private int numObjs;
    private int numTimeSteps;
    private ProblemSpec ps;
    private ArrayList<GObject> objects = new ArrayList<GObject>();

    public static void main(String[] args) {
        Vis v = new Vis();
    }

    public Vis(){
        super ("Visualiser");
        setDefaultCloseOperation (JFrame.EXIT_ON_CLOSE);
        prepareGUI();
        //mainFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }


    private void prepareGUI() {
        /*mainFrame = new JFrame("Visualiser");
        mainFrame.
        mainFrame.setLayout(new BorderLayout());
        mainFrame.addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent windowEvent) {
                System.exit(0);
            }
        });*/
        ps = new ProblemSpec();
        try {
            ps.loadProblem("output.txt");
        } catch (IOException e) {
            System.out.print(e);
        }
        if (!ps.problemLoaded()){
            System.out.print("\nOutput.txt not found in directory. Exiting... \n");
            System.exit(2);
        }
        timeStep = 1;
        numObjs = ps.numObjects;
        numTimeSteps = ps.numTimeSteps;

        interactionObject_ = null;
        color_             = null;

        // Create the GUI
        JPanel topLevel = new JPanel();
        topLevel.setLayout (new BorderLayout());
        getContentPane().add (topLevel);


        // Create the graphic canvas
        GWindow window = new GWindow();
        topLevel.add(window.getCanvas(), BorderLayout.CENTER);

        // Create scane with default viewport and world extent settings
        scene_ = new GScene (window);

        /*double w0[] = {-100, 100, 0.0};
        double w1[] = {-100, 100, 0.0};
        double w2[] = {0, 0, 0.0};
        scene_.setWorldExtent (w0, w1, w2);*/
        GWorldExtent e = scene_.getWorldExtent();

        // Create som graphic objects
        createGraphics(scene_);

        pack();
        setSize(new Dimension(800, 800));
        setVisible (true);

        // Install the selection interaction
        window.startInteraction (this);
        timeStep = 2;
        while (drawFlag){
            long startTime = System.currentTimeMillis();
            if (timeStep >= numTimeSteps){
                timeStep = 2;
            }
            for (int i = 0; i < numObjs; i++) {
                //ObjConfig configOld = ps.timesteps.get(timeStep-1).get(i);
                //TestObject o = (TestObject) scene_.find((int) configOld.xpos*100, (int) configOld.ypos*100);
                TestObject o = (TestObject) objects.get(i);
                ObjConfig config = ps.timesteps.get(timeStep).get(i);
                Point2D pos = new Point((int) config.xpos*3, (int) config.ypos*3);
                o.update(pos);
                o.refresh();
            }
            timeStep += 1;
            long endTime = System.currentTimeMillis();
            //System.out.println("Total time: " + (endTime - startTime) + "ms\n");
            try {
                long time = (40-(endTime - startTime));
                if (time > 0){
                    Thread.sleep(time);
                }
            } catch(InterruptedException ex) {
                Thread.currentThread().interrupt();
            }
        }
    }


    private void createGraphics (GScene scene)
    {
        int n = numObjs;
        for (int i = 0; i < n; i++) {
            Color color = Color.getHSBColor ((float) i / n, (float) 0.5, (float) 1.0);
            //Point2D pos = new Point((int) ((Math.random()*2-1)*250), (int) ((Math.random()*2-1)*250));
            ObjConfig config = ps.timesteps.get(timeStep).get(i);
            Point2D pos = new Point((int) config.xpos*3, (int) config.ypos*3);
            GObject object = new TestObject ("" + i, pos, color);
            objects.add(object);
            scene.add (object);
        }
    }



    /**
     * Handle button interactions.
     *
     * @param event  Event causing call to this method.
     */
    public void actionPerformed (ActionEvent event)
    {
        if (interactionObject_ == null)
            return;

        Object source = event.getSource();

        GObject parent = interactionObject_.getParent();

        if (source == frontButton_)
            parent.reposition (interactionObject_, parent.front());

        else if (source == backButton_)
            parent.reposition (interactionObject_, parent.back());

        else if (source == forwardButton_)
            parent.reposition (interactionObject_, parent.forward());

        else if (source == backwardButton_)
            parent.reposition(interactionObject_, parent.backward());

        frontButton_.setEnabled (!interactionObject_.isInFront());
        forwardButton_.setEnabled (!interactionObject_.isInFront());
        backButton_.setEnabled (!interactionObject_.isInBack());
        backwardButton_.setEnabled (!interactionObject_.isInBack());

        interactionObject_.refresh();
    }



    /**
     * Handle graphic events.
     *
     * @param event  Event type.
     * @param x,y    Cursor position.
     */
    public void event (GScene scene, int event, int x, int y)
    {

        // We care of button 1 clicks only
        if (event == GWindow.BUTTON1_UP) {

            GObject object = scene_.find (x, y);
            System.out.print("X=" + x + " - Y = " + y + "\n");

            if (object != null) {
                if (interactionObject_ != null)
                    interactionObject_.getStyle().setBackgroundColor (color_);

                interactionObject_ = object;
                color_ = interactionObject_.getStyle().getBackgroundColor();
                interactionObject_.getStyle().
                        setBackgroundColor (new Color (255, 255, 255));

                object.refresh();
            }
        }
    }



    /**
     * Defines the geometry and presentation for a sample
     * graphic object.
     */
    private class TestObject extends GObject
    {
        private Point2D pos;
        private GSegment  largeCircle_;
        private GSegment  smallCircle_;
        private GSegment  arm_;


        TestObject (String name, Point2D pos, Color color)
        {
            this.pos = pos;

            GStyle style = new GStyle();
            style.setBackgroundColor (color);
            style.setLineStyle (GStyle.LINESTYLE_INVISIBLE);
            setStyle (style);

            largeCircle_ = new GSegment();
            addSegment (largeCircle_);

        }


        public void draw()
        {
            // Center of viewport
            int x0     = (int) Math.round (getScene().getViewport().getCenterX());
            int y0     = (int) Math.round (getScene().getViewport().getCenterY());

            int width  = (int) Math.round (getScene().getViewport().getWidth());
            int height = (int) Math.round (getScene().getViewport().getHeight());

            int length  = Math.min (width, height) - 20;


            int[] largeCircleCoord = Geometry.createCircle (x0, y0,
                    (int) Math.round (length * 0.006));

            Matrix4x4 matrix = new Matrix4x4();

            matrix.setIdentity();
            matrix.translate (-x0, -y0);
            matrix.translate (+x0, +y0);
            matrix.translate (pos.getX(), pos.getY());

            matrix.transformXyPoints (largeCircleCoord);


            largeCircle_.setGeometry (largeCircleCoord);
        }

        public void update(Point2D pos){
            this.pos = pos;
            draw();
        }
    }


}
