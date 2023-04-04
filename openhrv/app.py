import sys
import logging
from PySide6.QtWidgets import QApplication
from openhrv.view import View
from openhrv.model import Model

logging.basicConfig(level=logging.DEBUG, handlers=[logging.StreamHandler()])
log = logging.getLogger()


class Application(QApplication):
    def __init__(self, sys_argv):
        super(Application, self).__init__(sys_argv)
        self._model = Model()
        self._view = View(self._model)


def main():
    app = Application(sys.argv)
    app._view.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
