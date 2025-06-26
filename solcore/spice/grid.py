"""Classes that generate a metalisation pattern for a solar cell.
"""


import numpy as np 
import pixie
from PIL import Image, ImageOps
import tempfile
from pathlib import Path


class GridPattern:
    """Representation of a metalisation pattern on the front surface of a solar cell.

    Instead of instantiating this class directly, subclass and implement the `draw` method to render a grid pattern.

    Discussion
    ----------
    You should use three grayscale pixel values to illustrate the solar cell metalization:
    - black (0.0), represents no metalisation
    - grey (0.5), represents grid fingers
    - white (1.0), represents the bus bar.
    """
    def draw(self) -> pixie.Image:
        raise NotImplementedError("The draw() method should be implemented by subclasses to draw specific grid patterns.")

    def save_as_image(self, path):
        image: pixie.Image = self.draw()
        image.write_file(path)  # This file is seems to be corrupt!
        img = Image.open(path)  # But, it can be opened by PIL and re-saved.

        # The shape of this array will be something like (300, 300, 4)
        # because pixie saves colours as RGBA. We need to conver this
        # to a gray scale image
        img = ImageOps.grayscale(img)
        img.save(path)
    
    def as_array(self) -> np.ndarray:
        # Write the image to a temporary directory,
        # load the image back using PIL and return 
        # the data as an array.
        with tempfile.TemporaryDirectory() as dirname:
            path = Path(dirname) / "grid.png"  # file name does matter
            self.save_as_image(path.as_posix())
            img = Image.open(path.as_posix())
            return np.asarray(img)

    @property
    def is_metal(self) -> np.ndarray:
        """Return a bool array where the pattern contains metal.
        """
        pattern = self.as_array()
        return np.where((pattern / pattern.max()) > 0.2, True, False)


class HGridPattern(GridPattern):
    """A classic H pattern concenrator solar cell pattern.
    """

    def __init__(self, bus_px_width, finger_px_widths, offset_px=5, nx=300, ny=300):
        self.bus_px_width = bus_px_width
        self.finger_px_widths = finger_px_widths
        self.offset_px = offset_px
        self.nx = nx
        self.ny = ny
    
    def draw(self) -> pixie.Image:

        bus_px_width = self.bus_px_width
        finger_px_widths = self.finger_px_widths
        offset_px = self.offset_px
        nx = self.nx
        ny = self.ny

        BUS_PAINT = pixie.Paint(pixie.SOLID_PAINT)
        BUS_PAINT.color = pixie.Color(1, 1, 1, 1)  # White

        FINGER_PAINT = pixie.Paint(pixie.SOLID_PAINT)
        FINGER_PAINT.color = pixie.Color(0.5, 0.5, 0.5, 1)  # Gray

        # Fill the image with black i.e. no metal. We are going 
        # to draw on top of this canvas using the BUS_PAINT and
        # the FINGER_PAINT paints.
        self.image = image = pixie.Image(nx, ny)
        BLACK = pixie.Color(0, 0, 0, 1)
        image.fill(BLACK)

        
        # NB Top-left corner is the (0, 0)
        ctx = image.new_context()
        ctx.fill_style = BUS_PAINT
        ctx.fill_rect(offset_px, offset_px, nx - 2 * offset_px, bus_px_width)
        ctx.fill_rect(offset_px , ny - offset_px - bus_px_width, nx - 2 * offset_px, bus_px_width)

        # The image now looks like this, with the bus bars drawn
        #
        #  ***************
        #  ***************
        #
        #
        #
        #
        #  ***************
        #  ***************
        #

        ctx = image.new_context()
        ctx.stroke_style = FINGER_PAINT
        f_origin_y = np.rint(bus_px_width + offset_px)
        f_length = np.rint(ny - 2 * offset_px - 2 * bus_px_width)
        n = len(finger_px_widths)
        w = nx - 2 * offset_px  # width of mask
        d = w / (2*n) # the half-spacing between fingers
        f_x = [offset_px + d] # location of first grid finger
        for idx in range(1, n):
            f_x.append(f_x[idx-1] + 2*d)
        
        # Round the x locations to the nearest integer, this avoid
        # the drawing tool kit from blending the colors
        f_x = np.rint(np.array(f_x))

        
        for f_origin_x, f_width in zip(f_x, finger_px_widths):

            ctx.stroke_segment(f_origin_x, f_origin_y, f_origin_x, f_origin_y + f_length)

        # The image now looks like this, with the n grid fingers drawn, here n = 3
        #
        #  ***************
        #  ***************
        #    |    |    |
        #    |    |    |
        #    |    |    |
        #    |    |    |
        #  ***************
        #  ***************
        #

        return image
    
    