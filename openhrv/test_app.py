from PySide6.QtCore import QObject, Signal, QTimer
from random import randrange
from utils import get_address_or_uuid


class MockBluetoothAddress:
    def toString(self):
        return "31:41:59:26:53:58"


class MockBluetoothUUID:
    def toString(self):
        return "{0229c9aa-7f96-421c-8862-d6c0d9f82121}"  # random from uuid.uuid4()


class MockSensor:
    def name(self):
        return "MockSensor"

    def address(self):
        return MockBluetoothAddress()

    def deviceUuid(self):
        return MockBluetoothUUID()


class MockSensorScanner(QObject):
    sensor_update = Signal(object)
    status_update = Signal(str)

    def scan(self):
        polar_sensors = [MockSensor()]
        self.sensor_update.emit(polar_sensors)
        self.status_update.emit(f"Found {len(polar_sensors)} sensor(s).")


class MockSensorClient(QObject):
    ibi_update = Signal(object)
    status_update = Signal(object)

    def __init__(self):
        super().__init__()

        self.timer = QTimer()
        self.timer.setInterval(1000)
        self.timer.timeout.connect(self.simulate_ibi)

    def connect_client(self, sensor):
        self.status_update.emit(
            f"Connecting to mock sensor at {get_address_or_uuid(sensor)}."
        )
        self.timer.start()

    def disconnect_client(self):
        self.status_update.emit("Disconnecting from mock sensor.")
        self.timer.stop()

    def simulate_ibi(self):
        self.ibi_update.emit(randrange(700, 1400))


def main():
    """Mock sensor classes.

    Mock classes need to replace their mocked counterparts in namespace before
    the latter are imported elsewhere:
    https://stackoverflow.com/questions/3765222/monkey-patch-python-class
    """
    from openhrv import sensor  # noqa

    sensor.SensorClient = MockSensorClient
    sensor.SensorScanner = MockSensorScanner

    from openhrv.app import main as mock_main  # noqa

    mock_main()


if __name__ == "__main__":
    main()
