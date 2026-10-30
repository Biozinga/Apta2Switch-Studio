"""The single selected-replica chart refreshes without mixing simulations."""
from copy import deepcopy
import json
import unittest

import pyarrow as pa

from tests.ui_helpers import FrenchAppTest
from tests.test_ffs_progress_charts import request, stage, state


SCRIPT = '''
import streamlit as st
from aptaswitch_studio.ffs_progress_charts import replica_rows, selected_replica_chart
states = st.session_state['states']
request = st.session_state['request']
from aptaswitch_studio.table_ui import render_readonly_table
render_readonly_table(replica_rows(states, request))
chart = selected_replica_chart(states[st.session_state['selected']], request)
if chart is not None:
    st.altair_chart(chart, width='stretch', key='selected_replica_interfaces')
'''


def rendered_rows(chart):
    return pa.ipc.open_stream(chart.proto.datasets[0].data.data).read_all().to_pylist()


class SelectedReplicaRefreshTests(unittest.TestCase):
    def app(self, states):
        app = FrenchAppTest.from_string(SCRIPT)
        app.session_state['states'] = states
        app.session_state['request'] = request()
        app.session_state['selected'] = 0
        app.run()
        self.assertFalse(app.exception)
        return app

    def test_flux_only_shows_source_counts_and_burn_table_without_empty_graph(self):
        app = self.app([state(0, phase='flux'), state(1, phase='equilibration')])
        self.assertFalse(app.get('vega_lite_chart'))
        table = app.get("table")[0].value
        self.assertEqual(table.loc[0, 'Avancement étape'], '9/200 franchissements')
        self.assertEqual(table.loc[1, 'Avancement étape'], '150000/200000 pas')
        self.assertEqual(table.loc[2, 'État'], 'En attente')

    def test_atomic_count_refresh_changes_dataset_and_keeps_fixed_axis(self):
        states = [state(stages=[stage(successes=1, failures=0)]), state(1)]
        app = self.app(deepcopy(states))
        self.assertEqual(len(app.get('vega_lite_chart')), 1)
        first = app.get('vega_lite_chart')[0]
        first_name = first.proto.datasets[0].name
        self.assertEqual(rendered_rows(first)[0]['successes'], 1)
        states[0]['stages'][0]['successes'] = 9
        app.session_state['states'] = deepcopy(states)
        app.run()
        self.assertFalse(app.exception)
        second = app.get('vega_lite_chart')[0]
        self.assertEqual(rendered_rows(second)[0]['successes'], 9)
        self.assertNotEqual(second.proto.datasets[0].name, first_name)
        self.assertEqual(json.loads(second.proto.spec)['encoding']['x']['scale']['domain'], [0, 1000])
        self.assertEqual(app.get("table")[0].value.loc[0, 'Succès interface'], '9')

    def test_selection_changes_only_chart_source_not_counts_in_other_replicas(self):
        states = [state(stages=[stage(successes=2)]), state(1, stages=[stage(successes=90)])]
        app = self.app(deepcopy(states))
        self.assertEqual(rendered_rows(app.get('vega_lite_chart')[0])[0]['successes'], 2)
        app.session_state['selected'] = 1
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(len(app.get('vega_lite_chart')), 1)
        self.assertEqual(rendered_rows(app.get('vega_lite_chart')[0])[0]['successes'], 90)
        self.assertEqual(app.get("table")[0].value['Succès interface'].tolist(), ['2', '90', '—', '—'])
        self.assertEqual(app.session_state['states'], states)


if __name__ == '__main__':
    unittest.main()
