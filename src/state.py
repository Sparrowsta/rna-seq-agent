
class PlanExecuteState(dict):
    """Plan-Execute状态管理类"""
    
    def __init__(self):
        super().__init__()
        self.update({
            "input": "",
            "plan": [],
            "past_steps": [],
            "response": "",
            "messages": [],
            "mode": "normal"
        })
