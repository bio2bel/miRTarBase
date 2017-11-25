# -*- coding: utf-8 -*-

import flask_admin
from flask import Flask
from flask_admin.contrib.sqla import ModelView

from bio2bel_mirtarbase.manager import Manager
from bio2bel_mirtarbase.models import *


def add_admin(app, session, **kwargs):
    admin = flask_admin.Admin(app, **kwargs)
    admin.add_view(ModelView(Interaction, session))
    admin.add_view(ModelView(Mirna, session))
    admin.add_view(ModelView(Target, session))
    admin.add_view(ModelView(Evidence, session))
    admin.add_view(ModelView(Species, session))
    return admin


def get_app(connection=None, url=None):
    app = Flask(__name__)
    manager = Manager(connection=connection)
    add_admin(app, manager.session, url=url)
    return app


if __name__ == '__main__':
    app = get_app()
    app.run(debug=True, host='0.0.0.0', port=5000)
