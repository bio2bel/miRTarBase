# -*- coding: utf-8 -*-

import flask_admin
from flask import Flask
from flask_admin.contrib.sqla import ModelView

from bio2bel_mirtarbase.manager import Manager
from bio2bel_mirtarbase.models import *


def add_admin_view(app, manager, url='/'):
    admin = flask_admin.Admin(app, url=url)
    admin.add_view(ModelView(Interaction, manager.session))
    admin.add_view(ModelView(Mirna, manager.session))
    admin.add_view(ModelView(Target, manager.session))
    admin.add_view(ModelView(Evidence, manager.session))
    return admin


def get_app(connection=None):
    app = Flask(__name__)
    manager = Manager(connection=connection)
    add_admin_view(app, manager)
    return app


if __name__ == '__main__':
    app = get_app()
    app.run(debug=True, host='0.0.0.0', port=5000)