# -*- coding: utf-8 -*-

import logging

from .manager import Manager

log = logging.getLogger(__name__)

main = Manager.get_cli()

if __name__ == '__main__':
    main()
