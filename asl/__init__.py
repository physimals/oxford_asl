# pkg-util style namespace package to allow asl package to be split over multiple locations
__path__ = __import__('pkgutil').extend_path(__path__, __name__)
