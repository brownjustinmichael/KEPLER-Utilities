import numpy as np
from numpy import ma
from matplotlib import scale as mscale
from matplotlib import transforms as mtransforms
from matplotlib.ticker import Formatter, FixedLocator
from matplotlib.ticker import (NullLocator, LogLocator, AutoLocator,
                               SymmetricalLogLocator, Locator, LogFormatter, NullFormatter)
from matplotlib.ticker import is_close_to_int, is_decade, decade_down, decade_up, nearest_long
from matplotlib import rcParams
import math

class ShiftLogLocator(Locator):
    """
    Determine the tick locations for log axes
    """

    def __init__(self, base=10.0, sign = 1.0, zero = 0.0, subs=[1.0], numdecs=4, numticks=15):
        """
        place ticks on the location= base**i*subs[j]
        """
        self.base(base)
        self.sign = sign
        self.zero = zero
        self.subs(subs)
        self.numticks = numticks
        self.numdecs = numdecs

    def base(self, base):
        """
        set the base of the log scaling (major tick every base**i, i integer)
        """
        self._base = base

    def subs(self, subs):
        """
        set the minor ticks the log scaling every base**i*subs[j]
        """
        if subs is None:
            self._subs = None  # autosub
        else:
            self._subs = np.asarray(subs)

    def __call__(self):
        'Return the locations of the ticks'
        vmin, vmax = self.axis.get_view_interval()
        return self.tick_values(vmin, vmax)

    def tick_values(self, vmin, vmax):
        b = self._base
        # dummy axis has no axes attribute
        if hasattr(self.axis, 'axes') and self.axis.axes.name == 'polar':
            vmax = math.ceil(math.log(vmax - self.zero) / math.log(b))
            decades = np.arange(vmax - self.numdecs, vmax)
            ticklocs = b ** decades

            return ticklocs
        if self.sign > 0.0:
            if vmin <= self.zero:
                if self.axis is not None:
                    vmin = self.axis.get_minpos()
                if vmax <= self.zero:
                    raise ValueError(
                        "Data has no positive values, and therefore can not be "
                        "log-scaled.")
        else:
            if vmax >= self.zero:
                if self.axis is not None:
                    vmin = self.axis.get_minpos()
                if vmin >= self.zero:
                    raise ValueError(
                        "Data has no positive values, and therefore can not be "
                        "log-scaled.")

        vmin = math.log((vmin - self.zero) * self.sign) / math.log(b)
        vmax = math.log((vmax - self.zero) * self.sign) / math.log(b)

        if vmax < vmin:
            vmin, vmax = vmax, vmin

        numdec = math.floor(vmax) - math.ceil(vmin)

        if self._subs is None:  # autosub
            if numdec > 10:
                subs = np.array([1.0])
            elif numdec > 6:
                subs = np.arange(2.0, b, 2.0)
            else:
                subs = np.arange(2.0, b)
        else:
            subs = self._subs

        stride = 1
        while numdec / stride + 1 > self.numticks:
            stride += 1

        decades = np.arange(math.floor(vmin) - stride,
                            math.ceil(vmax) + 2 * stride, stride)
        if hasattr(self, '_transform'):
            ticklocs = self._transform.inverted().transform(decades)
            if len(subs) > 1 or (len(subs == 1) and subs[0] != 1.0):
                ticklocs = np.ravel(np.outer(subs, ticklocs))
        else:
            if len(subs) > 1 or (len(subs == 1) and subs[0] != 1.0):
                ticklocs = []
                for decadeStart in b ** decades:
                    ticklocs.extend(subs * decadeStart)
            else:
                ticklocs = float (b) ** decades
        return self.raise_if_exceeds(np.asarray(ticklocs) * self.sign + self.zero)

    def view_limits(self, vmin, vmax):
        'Try to choose the view limits intelligently'
        b = self._base
        
        if vmax < vmin:
            vmin, vmax = vmax, vmin

        if self.axis.axes.name == 'polar':
            vmax = math.ceil(math.log(vmax - self.zero) / math.log(b))
            vmin = b ** (vmax - self.numdecs)
            return vmin, vmax

        minpos = self.axis.get_minpos()

        if (self.sign > 0.0) :
            if vmax <= self.zero or not np.isfinite(minpos):
                raise ValueError(
                    "Data has no positive values, and therefore can not be "
                    "log-scaled.")
        else:
            if vmin >= self.zero or not np.isfinite(minpos):
                raise ValueError(
                    "Data has no positive values, and therefore can not be "
                    "log-scaled.")
            
        if not is_decade((vmin - self.zero) * self.sign, self._base):
            if self.sign > 0.0:
                vmin = decade_down((vmin - self.zero) * self.sign, self._base) + self.zero
            else:
                vmin = decade_up((vmin - self.zero) * self.sign, self._base) * self.sign + self.zero
        if not is_decade((vmax - self.zero) * self.sign, self._base):
            if self.sign > 0.0:
                vmax = decade_up((vmax - self.zero) * self.sign, self._base) + self.zero
            else:
                vmax = decade_down((vmax - self.zero) * self.sign, self._base) * self.sign + self.zero
        if vmin == vmax:
            if (self.sign > 0.0):
                vmin = decade_down((vmin - self.zero) * self.sign, self._base) + self.zero
                vmax = decade_up((vmax - self.zero) * self.sign, self._base) + self.zero
            else:
                vmin = decade_up((vmin - self.zero) * self.sign, self._base) * self.sign + self.zero
                vmax = decade_down((vmax - self.zero) * self.sign, self._base) * self.sign + self.zero
            
        result = mtransforms.nonsingular(vmin, vmax)
        return result


class ShiftLogFormatterMathtext(LogFormatter):
    """
    Format values for log axis; using ``exponent = log_base(value)``
    """
    
    def __init__(self, base, sign, zero):
        LogFormatter.__init__ (self, base)
        self.sign = sign
        self.zero = zero
        # self.set_offset_string ("+%.f" % self.zero)
        
    def get_offset (self):
        usetex = rcParams['text.usetex']
        if self.zero < 0:
            sign_string = '-'
        elif self.zero > 0:
            sign_string = '+'
        else:
            return ''
        if usetex:
            return (r'$%s%.2e$' % (sign_string, abs (self.zero)))
        else:
            return ('$\mathdefault{%s%.2e}$' % (sign_string, abs (self.zero)))

    def __call__(self, x, pos=None):
        'Return the format for tick val *x* at position *pos*'
        b = self._base
        usetex = rcParams['text.usetex']

        # only label the decades
        if x == self.zero:
            if usetex:
                return '$0$'
            else:
                return '$\mathdefault{0}$'
                
        fx = math.log(abs(x - self.zero)) / math.log(b)
        is_decade = abs (round (fx) - fx) < 1.0e10

        sign_string = '-' if x < 0 else ''
        sign_string = '-' if self.sign < 0 else ''

        # use string formatting of the base if it is not an integer
        if b % 1 == 0.0:
            base = '%d' % b
        else:
            base = '%s' % b
        if not is_decade and self.labelOnlyBase:
            return ''
        elif not is_decade:
            if usetex:
                return (r'$%s%s^{%.2f}$') % \
                                            (sign_string, base, fx)
            else:
                return ('$\mathdefault{%s%s^{%.2f}}$') % \
                                            (sign_string, base, fx)
        else:
            if usetex:
                return (r'$%s%s^{%d}$') % (sign_string,
                                           base,
                                           nearest_long(fx))
            else:
                return (r'$\mathdefault{%s%s^{%d}}$') % (sign_string,
                                                         base,
                                                         nearest_long(fx))

class ShiftLogScale(mscale.ScaleBase):
    """
    Scales data in range -pi/2 to pi/2 (-90 to 90 degrees) using
    the system used to scale latitudes in a Mercator projection.

    The scale function:
      ln(tan(y) + sec(y))

    The inverse scale function:
      atan(sinh(y))

    Since the Mercator scale tends to infinity at +/- 90 degrees,
    there is user-defined threshold, above and below which nothing
    will be plotted.  This defaults to +/- 85 degrees.

    source:
    http://en.wikipedia.org/wiki/Mercator_projection
    """

    # The scale class must have a member ``name`` that defines the
    # string used to select the scale.  For example,
    # ``gca().set_yscale("mercator")`` would be used to select this
    # scale.
    name = 'shiftlog'


    def __init__(self, axis, **kwargs):
        """
        Any keyword arguments passed to ``set_xscale`` and
        ``set_yscale`` will be passed along to the scale's
        constructor.

        thresh: The degree above which to crop the data.
        """
        mscale.ScaleBase.__init__(self)
        base = kwargs.pop("base", 10.0)
        sign = kwargs.pop("sign", 1.0)
        if sign == 0.0:
            raise ValueError ("sign must not be zero")
        sign = sign / abs (sign)
        zero = kwargs.pop("zero", 0.0)
        if base <= 0.0:
            raise ValueError("base must be positive")
        self.base = base
        self.sign = sign
        self.zero = zero

    def get_transform(self):
        """
        Override this method to return a new instance that does the
        actual transformation of the data.

        The ShiftLogTransform class is defined below as a
        nested class of this one.
        """
        return self.ShiftLogTransform(self.base, self.sign, self.zero)

    def set_default_locators_and_formatters(self, axis):
        """
        Override to set up the locators and formatters to use with the
        scale.  This is only required if the scale requires custom
        locators and formatters.  Writing custom locators and
        formatters is rather outside the scope of this example, but
        there are many helpful examples in ``ticker.py``.

        In our case, the Mercator example uses a fixed locator from
        -90 to 90 degrees and a custom formatter class to put convert
        the radians to degrees and put a degree symbol after the
        value::
        """
        axis.set_major_locator(ShiftLogLocator(self.base, self.sign, self.zero))
        axis.set_major_formatter(ShiftLogFormatterMathtext(self.base, self.sign, self.zero))
        axis.set_minor_locator(ShiftLogLocator(self.base, self.sign, self.zero, []))
        axis.set_minor_formatter(NullFormatter())

    def limit_range_for_scale(self, vmin, vmax, minpos):
        """
        Override to limit the bounds of the axis to the domain of the
        transform.  In the case of Mercator, the bounds should be
        limited to the threshold that was passed in.  Unlike the
        autoscaling provided by the tick locators, this range limiting
        will always be adhered to, whether the axis range is set
        manually, determined automatically or changed through panning
        and zooming.
        """
        if self.sign > 0.0:
            return max(vmin, self.zero * (1.0 + 1.0e-15)), vmax
        else:
            return vmin, min (vmax, self.zero * (1.0 - 1.0e-15))
        # if self.sign > 0.0:
        #     return max(vmin, self.zero), vmax
        # else:
        #     return vmin, min (vmax, self.zero)
            
    class ShiftLogTransform(mtransforms.Transform):
        # There are two value members that must be defined.
        # ``input_dims`` and ``output_dims`` specify number of input
        # dimensions and output dimensions to the transformation.
        # These are used by the transformation framework to do some
        # error checking and prevent incompatible transformations from
        # being connected together.  When defining transforms for a
        # scale, which are, by definition, separable and have only one
        # dimension, these members should always be set to 1.
        input_dims = 1
        output_dims = 1
        is_separable = True

        def __init__(self, base, sign, zero):
            mtransforms.Transform.__init__(self)
            self.base = base
            self.sign = sign
            self.zero = zero

        def transform_non_affine(self, a):
            """
            This transform takes an Nx1 ``numpy`` array and returns a
            transformed copy.  Since the range of the Mercator scale
            is limited by the user-specified threshold, the input
            array must be masked to contain only valid values.
            ``matplotlib`` will handle masked arrays and remove the
            out-of-range data from the plot.  Importantly, the
            ``transform`` method *must* return an array that is the
            same shape as the input array, since these values need to
            remain synchronized with values in the other dimension.
            """
            masked = ma.masked_where(((a - self.zero) * self.sign <= 0.0), (a - self.zero) * self.sign)
            if masked.mask.any():
                return ma.log10(masked) / np.log10 (self.base)
            else:
                return np.log10((a - self.zero) * self.sign) / np.log10 (self.base)

        def inverted(self):
            """
            Override this method so matplotlib knows how to get the
            inverse transform for this transform.
            """
            return ShiftLogScale.InvertedShiftLogTransform(self.base, self.sign, self.zero)

    class InvertedShiftLogTransform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True

        def __init__(self, base, sign, zero):
            mtransforms.Transform.__init__(self)
            self.base = base
            self.sign = sign
            self.zero = zero

        def transform_non_affine(self, a):
            return np.power(self.base, a) * self.sign + self.zero

        def inverted(self):
            return ShiftLogScale.ShiftLogTransform(self.base, self.sign, self.zero)

# Now that the Scale class has been defined, it must be registered so
# that ``matplotlib`` can find it.
mscale.register_scale(ShiftLogScale)